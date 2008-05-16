/* MODULE:      MMG
 *
 * FILE:	MMG_make_dot
 *
 * DESCRIPTION:	Outputs a DOT file with the network in it.
 *
 * VERSION:     1.2.1
 *
 * AUTHOR:      Josselin Noirel <j.noirel@sheffield.ac.uk>
 *
 * REFERENCES:  Sanguinetti et al., Bioinformatics (2008).
 *              MMG: a probabilistic tool to identify submodules of metabolic pathways.
 *
 * NOTE:        Based on Guido Sanguinetti's code available from
 *	        http://www.dcs.shef.ac.uk/~guido/.
 *
 * LICENCE:     GPL-3
 *
 * NOTES:	* One tricky thing is the fact that the numbering
 *		in the data file is 1 .. N, whereas here it is the
 *		C convention, which is used 0 .. N-1.  Macros could
 *		have been devised READ_FROM_R() and PRINT_TO_R().
 *		As for now, the instances where such manipulations are
 *		needed are flagged with 'C-numbering'.
 *
 *		Basically, as soon as possible everything is done
 *		using the C convention.  Turning back to the R
 *		convention is done only at the last minute at the time
 *		of outputting something.
 *
 *		* Also, the matrices are flattened by R but in a
 *              surprising way.  {{a, b}, {c, d}} becomes {a, c, b, d}.
 *
 *		* It is to be noted as well that the adjacency lists
 *              describe incoming edges.  The neighbours of i
 *              influence i!  And not: The neighbours of i are
 *              influenced by i.  This is because the class of i is
 *              influenced by its neighbours and not the opposite.
 *
 * ACKNOWLEDGMENTS: The author expresses his gratitude to Robert
 *		Sochon (University of Sheffield, the UK) and William
 *		Venables (CSIRO, Australia).
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>

#define DEBUG 0 /* Only small networks! */
#define R     1 /* Comply with Bioconductor's requirements */

/* DECLARATIONS */

void	       MMG_compute (int *, int *, double *, int *, double *, double *,
			    double *, int *, int *, int *, double *, double *);

void	 MMG_Gibbs_sampler (unsigned, size_t, size_t, int *, double *, int
			    *, double *, int (*)[3], int *,
			    double *, double *, double, double, double *, double *,
			    double, double);

double	      sum_log_rand (unsigned);

void   MMG_class_posterior (size_t, size_t, int *, double *, int *,
			    double *, int (*)[3], double, double,
			    double, double, double [3]);

/* FUNCTIONS */

/* MMG_compute()
 *
 *	The matlab equivalent is MMGmasterSampler.m.
 *
 *	Function intended to be .C-called from within R.  Some
 *	variables are modified in place: - samples, - lup, and - ldown.
 *	This function does the following, given the parameters and the
 *	data:
 *
 *	1. Allocate classes[].
 *
 *	2. Initialise classes[] depending on the data.
 *
 *	3. Initialise lambda_up, lambda_down.
 *
 *	4a. Burn-in;
 *	4b. Simulation.
 *
 *	At stage 4a-b, one calls MMG_Gibbs_sampler.
 *
 * NB:	There is an important difference to take into consideration
 *	regarding the variables nbs[], wts[], and samples[].  Whilst
 *	nbs[] and wts[] are matrices in R passed as vectors (and
 *	flattened according to R rules), samples[] was declared as a
 *	vector from the very start (even though it is conceptually a
 *	matrix n x 3).
 *
 *	Therefore the conceptual elements nbs[i, j] and wts[i, j] are
 *	accessed using nbs[i + j * n] and wts[i + j * n], and
 *	samples[i, j] using samples[i * 3 + j].  (n is the number of
 *	nodes).
 */

void MMG_compute (int	 * dims,    /* dims[0] x dims[1] */
		  int	 * nbnums,  /* nbnums[i]	 i < dims[0] */
		  double * msvals,  /* msvals[i]	 i < dims[0]  NB: Can be NaN */
		  int	 * nbs,	    /* nbs[i, j]	 i < dims[0], j < nbnums[i]  */
		  double * wts,	    /* wts[i, j]	 i < dims[0], j < nbnums[i]  */
		  double * sigma,   /* Parameters of the model (fixed) */
		  double * alpha,   /* : */
		  int	 * burn_in, /* Burn-in (if one considers this useful) */
		  int	 * steps,   /* Number of steps of simulation */
		  int	 * samples, /* samples[i, 3]  i < dims[0] */
		  double * lup,	    /* lup[i]	      i < steps	  */
		  double * ldown)   /* ldown[i]	      i < steps	  */
{
	size_t n, m;		       /* Number of nodes, maximum number of neighbours */
	double mean, mean2, stddev;    /* Average of the data */
	size_t n_up, n_down;	       /* Number of up/down */
	double lambda_up, lambda_down; /* Cf. article */
	double s_up, s_down;           /* Sum of {up, down}-data */
	double min, max;	       /* min, max of the data */
	int (*classes)[3];             /* Classes at time t */

	/* 0. Initialisations */

	n = (size_t)dims[0];
	m = (size_t)dims[1];
	classes = NULL;

	/* 1. Allocate the arrays */

#if R
	/* classes = Calloc(n, int[3]); B***y macro! It needs no casting... */
	classes = R_chk_calloc(n, sizeof(int[3]));
#else
	classes = malloc(n * sizeof(int[3]));
#endif
	if (classes == NULL) {
		error("Could not allocate memory for classes.");
		goto RETURN;
	}

	for (size_t i = 0; i < n; i++) {
		classes[i][0] = classes[i][1] = classes[i][2] = 0;
		samples[i * 3] = samples[i * 3 + 1] = samples[i * 3 + 2] = 0;
	}

	/* 2. Compute the average, the standard deviation, min, max */
	      /* NA considered as zero */
              /* We ensure that min <= 0.0 <= max (***) */

	min = max = 0.0; /* <-- (***) */
	mean = mean2 = 0.0;

	for (size_t i = 0; i < n; i++) {
		double x = msvals[i];
		if (! ISNA(x)) {
			mean  += x;
			mean2 += x * x;
			if (x < min) {
				min = x;
			}
			if (x > max) {
				max = x;
			}
		}
	}

	mean  /= (double)n; /* n takes into account NaN */
	mean2 /= (double)n; /* Ditto */
	stddev = sqrt(mean2 - mean * mean);

	    /* It would be good to ensure that min < 0 and max > 0 */

	assert(min <= 0.0);
	assert(max >= 0.0);

	/* 3. Initialise the classes and the lambda* */
	      /* NA considered as zero */

	s_up = s_down = 0.0;
	n_up = n_down = 0;
	lambda_up = lambda_down = 0.0;

	for (size_t i = 0; i < n; i++) {
		double x = msvals[i];

		if	(! ISNA(x) && x < mean - 0.5 * stddev) {
			classes[i][0] = 1;
			s_down += x;
			n_down++;
		}
		else if (! ISNA(x) && x > mean + 0.5 * stddev) {
			classes[i][2] = 1;
			s_up += x;
			n_up++;
		}
		else { /* NA or unchanged */
			classes[i][1] = 1;
		}
	}

	if (n_up > 0) {
		s_up /= (double)n_up;
		lambda_up = +1.0 / s_up;
	}
	else {
		lambda_up = max + 2.0;
	}

	if (n_down > 0) {
		s_down /= (double)n_down;
		lambda_down = -1.0 / s_down;
	}
	else {
		lambda_down = 2.0 - min;
	}

	assert(lambda_up   > 0.0);    /* The lambdas must be > 0 */
	assert(lambda_down > 0.0);

	/* 4. Run the Gibbs sampler */
	    /* Burn-in first */

#if R
	GetRNGstate();
#endif

	for (unsigned t = 0U; t < (unsigned)burn_in[0]; t++) {
		MMG_Gibbs_sampler(t, n, m, nbnums, msvals, nbs, wts,	   /* Data */
				  classes, NULL, NULL, NULL,		   /* Arrays */
				  sigma[0], alpha[0], &lambda_up, &lambda_down, min, max); /* Parameters */
	}

	    /* Collect data */

	for (unsigned t = 0U; t < (unsigned)steps[0]; t++) {
		MMG_Gibbs_sampler(t, n, m, nbnums, msvals, nbs, wts,	   /* Data */
				  classes, samples, lup, ldown,		   /* Arrays */
				  sigma[0], alpha[0], &lambda_up, &lambda_down, min, max); /* Parameters */
	}

#if R
	PutRNGstate();
#endif

RETURN:
#if R
#	define FREE(p) Free(p)
#else
#	define FREE(p) free(p)
#endif

	if (classes != NULL)
		FREE(classes);

#undef FREE

	return;
}

/* FUNCTION: MMG_Gibbs_sampler()
 *
 *	The matlab equivalent is MMGGibbsSampler.m.
 *
 *	Computes the posterior probability f[] using
 *	MMG_class_posterior() and update the classes, for each node.
 *	Also, lambda_up and lambda_down are updated (to that purpose,
 *	they are passed as pointers).
 *
 * NB:	The main difference there is between this implementation and
 *	the Matlab implementation is that some of the code used to
 *	update lambdaUp and lambdaDown are moved into the loop
 *	responsible for the class update.
 */

#if R
#	define RAND1() unif_rand()
#else
#	define RAND1() ((double)rand() / (double)RAND_MAX) /* X = Unif([0,1]) */
#endif

void MMG_Gibbs_sampler (unsigned t,	      /* Step */
			size_t n, size_t m,   /* Dimensions */
			int    * nbnums,      /* Number of neighbours */
			double * msvals,      /* MS measurements */
			int    * nbs,	      /* Neighbours' indices */
			double * wts,	      /* Weights */
			int (*classes)[3],    /* Classes at time t */
			int *samples,	      /* Samples	       (can be NULL) */
			double * lup,	      /* Record of lambda_up   (can be NULL) */
			double * ldown,	      /* Record of lambda_down (can be NULL) */
			double sigma,	      /* Parameter */
			double alpha,	      /* Ditto */
			double * lambda_up,   /* Ditto NB: pointer because it has to be updated */
			double * lambda_down, /* Ditto NB: pointer because it has to be updated */
			double	 min,	      /* Range of the MS distribution */
			double	 max)
{
	unsigned n_up, n_down;
	double	 s_up, s_down; /* Sums [sum(dataUp/Down)] */

	n_up = n_down = 0U;  /* Number of ups and downs */
	s_up = s_down = 0.0; /* Sum of the corresponding data */

	for (size_t i = 0; i < n; i++) {
		double f[3];
		double x;

		/* This fills up f[] with the posterior probabilities */
		MMG_class_posterior(i, n, nbnums, msvals, nbs, wts, classes,
				    sigma, alpha, *lambda_up, *lambda_down, f);

		assert(f[0] >= 0.0 && f[0] <= 1.0);
		assert(f[1] >= 0.0 && f[1] <= 1.0);
		assert(f[2] >= 0.0 && f[2] <= 1.0);

		x = RAND1();

		if	(x < f[0]) {
			classes[i][0] = 1;
			classes[i][1] = 0;
			classes[i][2] = 0;

			n_down++;
			s_down += (ISNA(msvals[i]) ? -1.0 / *lambda_down : msvals[i]);

			if (samples != NULL)
				samples[3 * i]++;
		}
		else if (x < f[0] + f[1]) {
			classes[i][0] = 0;
			classes[i][1] = 1;
			classes[i][2] = 0;

			if (samples != NULL)
				samples[3 * i + 1]++;
		}
		else {
			classes[i][0] = 0;
			classes[i][1] = 0;
			classes[i][2] = 1;
			
			n_up++;
			s_up += (ISNA(msvals[i]) ? 1.0 / *lambda_up : msvals[i]);

			if (samples != NULL)
				samples[3 * i + 2]++;
		}
	}

	if (n_up > 0U) {
		*lambda_up = -sum_log_rand(n_up) / s_up;
	}
	else {
		*lambda_up = 1.0 / (max + sigma);
	}

	if (n_down > 0U) {
		*lambda_down = sum_log_rand(n_down) / s_down;
	}
	else {
		*lambda_down = 1.0 / (sigma - min);
	}

	assert(*lambda_up   > 0.0);
	assert(*lambda_down > 0.0);

	if (lup != NULL)
		lup[t] = *lambda_up;

	if (ldown != NULL)
		ldown[t] = *lambda_down;

	return;
}

/* sum_log_rand( N )
 *
 *	Draw a random number X from a Gamma distribution
 *
 *		X = \sum_{i = 1 .. N} \log U ,
 *
 *	where U is a random variable uniformly distributed on [0, 1].
 */

double sum_log_rand (unsigned n)
{
	double x = 0.0;

	for (unsigned i = 0U; i < n; i++) {
		x += log(RAND1());
	}

	return x;
}

#undef RAND1

/* MMG_class_posterior()
 * 
 *	Matlab equivalent MMGclassPosterior.m.
 *
 *      The function computes the posterior probabilities for a node I
 *	to be down-regulated, up-regulated, or unchanged.  It does:
 *
 *      1. It constructs a vector T = {t1, t2, t3} where
 *
 *         - t1 is the sum of the weights of the connections coming
 *         from down-regulated nodes,
 *
 *         - t2 is the sum of the weights of the connections coming
 *         from unchanged nodes,
 *
 *         - t3 is the sum of the weights of the connections coming
 *         from up-regulated nodes,
 *
 *	2. T is normalised: T := T / Î£T.
 *
 *	3. F = {f1, f2, f3} where
 *
 *         If data > 0 (up-regulated)
 *
 *         	f1 = 0
 *         	f2 = Gaussian(0, sigma)(data)
 *         	f3 = Exponential(lambda_up)(data)
 *
 *         If data < 0 (down-regulated)
 *
 *         	f1 = Exponential(lambda_down)(-data)
 *         	f2 = Gaussian(0, sigma)(data)
 *         	f3 = 0
 *
 *         If data == 0 || undefined
 *
 *              f1 = f2 = f3 = 1
 *              
 *	4. F := F . T = {f1 t1, f2 t2, f3 t3}
 *
 *      5. F is "normalised": F := F / Î£F.
 *
 */

void MMG_class_posterior (size_t ind, size_t n,
			  int	 * nbnums,
			  double * msvals,
			  int	 * nbs,
			  double * wts,
			  int (*classes)[3],
			  double sigma,
			  double alpha,
			  double lambda_up,
			  double lambda_down,
			  double f[3])
{
	double t[3], st, sf;   /* Intermediary */
	double d, dd, ss, sq;  /* Data */

	d  = msvals[ind];
	dd = msvals[ind] * msvals[ind];
	ss = sigma * sigma;
	sq = sqrt(2.0 * M_PI);

	t[0] = t[1] = t[2] = alpha;

	/* Loop over ind's neighbours */

	for (size_t i = 0; i < nbnums[ind]; i++) {
		size_t nb = (size_t)nbs[ind + i * n] - 1; /* C-numbering! */
		double wt =	    wts[ind + i * n];

		t[0] += wt * classes[nb][0];
		t[1] += wt * classes[nb][1];
		t[2] += wt * classes[nb][2];
	}

	st = t[0] + t[1] + t[2];
	t[0] /= st;
	t[1] /= st;
	t[2] /= st;

	if (d > 0) {
		f[0] = 0.0;
		f[1] = exp(-0.5 * dd / ss) / (sq * sigma);
		f[2] = lambda_up * exp(-lambda_up * d);
	}
	else if (d < 0) {
		f[0] = lambda_down * exp(lambda_down * d); /* No minus sign, d < 0 */
		f[1] = exp(-0.5 * dd / ss) / (sq * sigma);
		f[2] = 0.0;
	}
	else { /* NA */
		f[0] = 1.0;
		f[1] = 1.0;
		f[2] = 1.0;
	}

	f[0] *= t[0];
	f[1] *= t[1];
	f[2] *= t[2];

	sf = f[0] + f[1] + f[2];
	f[0] /= sf;
	f[1] /= sf;
	f[2] /= sf;

	return;
}
