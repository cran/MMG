/* MODULE:      MMG
 *
 * FILE:	MMG_make_dot
 *
 * DESCRIPTION:	Outputs a DOT file with the network in it.
 *
 * VERSION:     1.2.2
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
 * TODOs:	* Allow to set a limit for the weight.  Add two arguments.
 *		wmax and wmin (caution, the meanings with respect to
 *		MSKEGG are reverted).
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
#define R     1

/* DECLARATIONS */

void MMG_make_dot (int    *,
		   int    *,
		   int    *,
		   char  **,
		   int    *,
		   int    *,
		   int    *,
		   int    *,
		   double *);

double MMG_xlogx (double);

/* FUNCTIONS */

void MMG_make_dot (int    * dims,
		   int    * select,    /* List of integers */
		   int    * selectn,   /* Length of select */
		   char  ** out,       /* Output file */
		   int    * nbnums,
		   int    * nbs,
		   int    * rm_loops,  /* Should we avoid loops */
		   int    * type,      /* Directed (1) / Undirected (2)*/
		   double * samples)
{
#define NEIGHBOUR(i, j) (nbs[(i) + n * (j)] - 1) /* C-numbering */

	size_t n, m, p; /* Nodes, Max neighbours, Length selection */
	bool directed;  /* Directed ->, or Undirected -- graph */
	bool rmloops;   /* rm_loops but a more idiomatic type */
	size_t  * undirected_nbnums; /* nbnums */
	size_t ** undirected_nbs;    /* nbs */
	FILE    * fh;
	bool    * in_list;

	/* 0. Initialisations */

	n = (size_t)dims[0];
	m = (size_t)dims[1];
	p = (size_t)selectn[0];
	undirected_nbnums = NULL;
	undirected_nbs = NULL;
	in_list = NULL;
	directed = true;
	fh = NULL;

	if (rm_loops[0])
		rmloops = true;
	else
		rmloops = false;

	if (type[0] == 1)
		/* Directed */
		directed = true;
	else
		/* Undirected */
		directed = false;

#if R
	in_list = Calloc(n, bool);
#else
	in_list = malloc(n * sizeof(bool));
#endif

	if (in_list == NULL) {
		error("Could not allocate memory for in_list.\n");
		goto RETURN;
	}

	for (size_t i = 0U; i < n; i++)
		in_list[i] = false;

	for (size_t i = 0U; i < p; i++) {
		size_t sel = select[i] - 1U; /* C-numbering */
		in_list[sel] = true;
	}

	if (! directed) {
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
				error("Could not allocate undirected_nbs[].\n"); /* TODO */
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
	} /* End of if (! directed) */

	fh = fopen(out[0], "w");
	if (fh == NULL) {
		error("Could not open file '%s'\n", out[0]);
		goto RETURN;
	}

	if (directed)
		fprintf(fh, "digraph network\n{\n");
	else
		fprintf(fh, "graph network\n{\n");

	if (directed) {
		for (size_t i = 0U; i < n; i++) {
			if (! in_list[i])
				continue;

			for (size_t j = 0U; j < nbnums[i]; j++) {
				size_t nb = NEIGHBOUR(i, j);

				if (! in_list[nb])
					continue;

				if (! rmloops || nb != i)
					/* C-numbering */
					fprintf(fh, "\t%zu -> %zu;\n", nb + 1U, i + 1U); /* j -> i NOT i -> j */
			}
		}
	}
	else {
		for (size_t i = 0U; i < n; i++) {
			if (! in_list[i])
				continue;

			for (size_t j = 0U; j < undirected_nbnums[i]; j++) {
				size_t nb = undirected_nbs[i][j];

				if (! in_list[nb])
					continue;

				if (nb < i || (! rmloops && nb == i))
					/* C-numbering */
					fprintf(fh, "\t%zu -- %zu;\n", nb + 1U, i + 1U); /* j -> i */
			}
		}
	}

	/* Add a little bit of colour */

#define SAMPLE(i, j) (samples[(i) + n * (j)]) /* TODO: CHECK THIS AND LOCALISE */
#define MAKEHEX(x) MAKEHEY((int)(x * 256.0))
#define MAKEHEY(y) ((y) == 256 ? (y) - 1 : (y))

	for (size_t i = 0U; i < p; i++) {
		double pdown, punch, pup;
		size_t nd;

		nd = select[i] - 1U; /* C-numbering */

		pdown = SAMPLE(nd, 0);
		punch = SAMPLE(nd, 1);
		pup   = SAMPLE(nd, 2);

		fprintf(fh, "\t%zu [color=\"#%02x%02x%02x\", style=filled];\n", nd + 1U,
			MAKEHEX(pup), MAKEHEX(punch), MAKEHEX(pdown)); /* C-numbering */
	}

#undef SAMPLE
#undef MAKEHEX
#undef MAKEHEY

	fprintf(fh, "}\n");

	fclose(fh);
	fh = NULL;

RETURN:
	if (fh != NULL)
		fclose(fh);

#if R
#	define FREE(p) Free(p)
#else
#	define FREE(p) free(p)
#endif

	if (in_list != NULL)
		FREE(in_list);

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
