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
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <R.h>

#define DEBUG 0
#define R     0

/* DECLARATIONS */

void MMG_cut_graph (int    *,
		    int    *,
		    double *,
		    int    *,
		    int    *, /* Not used yet */
		    double *,
		    int    *,
		    int    *,
		    double *,
		    int    *,
		    int    *);

int MMG_find_connected_components (size_t,
				   int    *,
				   size_t *,
				   size_t **,
				   double *,
				   double,
				   int,
				   bool (*)(double, double, double, double, int));

void MMG_visit_node (size_t, size_t, int,
		     int    *,
		     size_t *,
		     size_t **,
		     double *,
		     double,
		     int,
		     bool  (*)(double, double, double, double, int));

bool MMG_method_likeliest (double, double, double, double, int);

bool MMG_method_threshold (double, double, double, double, int);

double MMG_xlogx (double);

bool MMG_method_entropy (double, double, double, double, int);

/* FUNCTIONS */

void MMG_cut_graph (int    * dims,
		    int    * nbnums,
		    double * msvals,
		    int    * nbs,
		    int    * wts,
		    double * samples,
		    int    * method,
		    int    * select,
		    double * threshold,
		    int    * category,
		    int    * ncateg)
{
#define NEIGHBOUR(i, j) (nbs[(i) + n * (j)] - 1) /* C-numbering */

	size_t n, m;
	size_t * undirected_nbnums; /* nbnums */
	size_t **undirected_nbs;    /* nbs */

	/* 0. Initialisations */

	n = (size_t)dims[0];
	m = (size_t)dims[1];
	undirected_nbnums = NULL;
	undirected_nbs = NULL;

	/* I. Build the undirected graph */

	/* 1. Allocation */

#if R
	undirected_nbnums = Calloc(n, size_t);
#else
	undirected_nbnums = malloc(n * sizeof(size_t));
#endif

	if (undirected_nbnums == NULL) {
		error("Could not allocate undirected_nbnums!\n");
		goto RETURN;
	}

#if R
	undirected_nbs = Calloc(n, size_t *);
#else
	undirected_nbs = malloc(n * sizeof(size_t *));
#endif

	if (undirected_nbs == NULL) {
		error("Could not allocate undirected_nbs!\n");
		goto RETURN;
	}

	for (size_t i = 0U; i < n; i++)
		undirected_nbs[i] = NULL;

	/* 2. Compute the number of neighbours */
	    /* First, initialise to the number of outward neighbours */

	for (size_t i = 0U; i < n; i++) {
		undirected_nbnums[i] = (size_t)nbnums[i];
	}

	    /* Second, add up the number of inward neighbours (that
	     * are not already counted in the outward ones) */

	for (size_t i = 0U; i < n; i++) {
		for (size_t j = 0U; j < (size_t)nbnums[i]; j++) {
			/* i <- nb */
			size_t nb = NEIGHBOUR(i, j);

			/* Was already i counted in the direction nb <- i? */
			bool found = false;
			for (size_t k = 0U; k < (size_t)nbnums[nb]; k++) {
				size_t nbnb = NEIGHBOUR(nb, k);

				if (nbnb == i) {
					found = true;
					break;
				}
			}

			if (! found)
				undirected_nbnums[nb]++;
		}
	}

	/* 3. Adjacency lists */
	    /* Allocate */

	for (size_t i = 0U; i < n; i++) {
#if R
		undirected_nbs[i] = Calloc(undirected_nbnums[i], size_t);
#else
		undirected_nbs[i] = malloc(undirected_nbnums[i] * sizeof(size_t));
#endif

		if (undirected_nbs[i] == NULL) {
			error("Could not allocate undirected_nbs[i]");
			goto RETURN;
		}

		undirected_nbnums[i] = 0U; /* Reinitialise */
	}

	    /* Basically, the same code as (2) but see the arrow <-- */

	for (size_t i = 0U; i < n; i++) {
		for (size_t j = 0U; j < (size_t)nbnums[i]; j++) {
			/* i <- nb */
			size_t nb = NEIGHBOUR(i, j);

			undirected_nbs[i][ undirected_nbnums[i]++ ] = nb;
		}
	}

	for (size_t i = 0U; i < n; i++) {
		for (size_t j = 0U; j < (size_t)nbnums[i]; j++) {
			/* i <- nb */
			size_t nb = NEIGHBOUR(i, j);

			/* Was already i counted in the direction nb <- i? */
			bool found = false;
			for (size_t k = 0U; k < (size_t)nbnums[nb]; k++) {
				size_t nbnb = NEIGHBOUR(nb, k);

				if (nbnb == i) {
					found = true;
					break;
				}
			}

			if (! found)
				undirected_nbs[nb][ undirected_nbnums[nb]++ ] = i; /* <-- */
		}
	}

	/* II. Find the connected components */

	ncateg[0] = MMG_find_connected_components(n, category,
						  undirected_nbnums,
						  undirected_nbs,
						  samples,
						  threshold[0],
						  select[0],
						  (method[0] == 1 ?  MMG_method_likeliest :
						   (method[0] == 2 ? MMG_method_threshold :
						                     MMG_method_entropy)));

RETURN:
#if R
#	define FREE(p) Free(p)
#else
#	define FREE(p) free(p)
#endif

	if (undirected_nbnums != NULL)
		FREE(undirected_nbnums);

	if (undirected_nbs != NULL) {
		for (size_t i = 0U; i < n; i++)
			if (undirected_nbs[i] != NULL)
				FREE(undirected_nbs[i]);
		FREE(undirected_nbs);
	}

#undef FREE

	return;

#undef NEIGHBOUR
}

 /* MMG_find_connected_components()
  *
  *	If there are N connected components, it will
  *	return N + 1.  The categories will therefore
  *	be 0 (uncertain), 1, 2 . . . N.
  */

#define SAMPLE(i, j) (samples[(i) + n * (j)]) /* TODO: CHECK THIS AND LOCALISE */

int MMG_find_connected_components (size_t n,
				   int    * category,
				   size_t * nbnums,
				   size_t **nbs,
				   double * samples,
				   double   threshold,
				   int      select,
				   bool  (* funct)(double, double, double, double, int))
{
	int cat = 1;

	for (size_t i = 0U; i < n; i++) {
		if (category[i] == 0 && funct(SAMPLE(i, 0), SAMPLE(i, 1), SAMPLE(i, 2), threshold, select)) {
			MMG_visit_node(i, n, cat, category, nbnums, nbs, samples, threshold, select, funct);
			cat++;
		}
	}

	return cat;
}

void MMG_visit_node (size_t ind,
		     size_t n,
		     int cat,
		     int    * category,
		     size_t * nbnums,
		     size_t **nbs,
		     double * samples,
		     double   threshold,
		     int      select,
		     bool  (* funct)(double, double, double, double, int))
{
	category[ind] = cat;

	for (size_t j = 0U; j < nbnums[ind]; j++) {
		size_t nb = nbs[ind][j];

		if (category[nb] == 0 && funct(SAMPLE(nb, 0), SAMPLE(nb, 1), SAMPLE(nb, 2), threshold, select))
			MMG_visit_node(nb, n, cat, category, nbnums, nbs, samples, threshold, select, funct);
	}

	return;
}

#undef SAMPLE

bool MMG_method_likeliest (double pdown, double punch, double pup, /* Posteriors */
			   double threshold, int select) /* threshold > 0 */
{
	bool res = false;

	if      (select == 1 && pdown > punch && pdown > pup   && pdown >= threshold)
		res = true;
	else if (select == 2 && punch > pdown && punch > pup   && punch >= threshold)
		res = true;
	else if (select == 3 && pup   > pdown && pup   > punch && pup   >= threshold)
		res = true;

	return res;
}

bool MMG_method_threshold (double pdown, double punch, double pup, /* Posteriors */
			   double threshold, int select) /* threshold > 0 */
{
	bool res = false;

	if      (select == 1 && pdown - punch >= threshold && pdown - pup   >= threshold)
		res = true;
	else if (select == 2 && punch - pdown >= threshold && punch - pup   >= threshold)
		res = true;
	else if (select == 3 && pup   - pdown >= threshold && pup   - punch >= threshold)
		res = true;

	return res;
}

double MMG_xlogx (double x) /* x >= 0.0 */
{
	return (x < DBL_EPSILON ? 0.0 : x * log(x));
}

bool MMG_method_entropy (double pdown, double punch, double pup, /* Posteriors */
			   double threshold, int select) /* threshold > 0 */
{
	bool res = false;

	double e = -(MMG_xlogx(pdown) + MMG_xlogx(punch) + MMG_xlogx(pup));

	if      (select == 1 && pdown > punch && pdown > pup   && e <= threshold)
		res = true;
	else if (select == 2 && punch > pdown && punch > pup   && e <= threshold)
		res = true;
	else if (select == 3 && pup   > pdown && pup   > punch && e <= threshold)
		res = true;

	return res;
}
