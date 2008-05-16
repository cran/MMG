# THIS FILE IS ENCODED IN -*- utf8 -*-

# MODULE:      MMG
#
# FILE:	MMG_make_dot
#
# DESCRIPTION:	Outputs a DOT file with the network in it.
#
# VERSION:     1.2.1
#
# AUTHOR:      Josselin Noirel <j.noirel@sheffield.ac.uk>
#
# REFERENCES:  Sanguinetti et al., Bioinformatics (2008).
#              MMG: a probabilistic tool to identify submodules of metabolic pathways.
#
# NOTE:        Based on Guido Sanguinetti's code available from
#	        http://www.dcs.shef.ac.uk/~guido/.
#
# LICENCE:     GPL-3
#
# NOTES:	* One tricky thing is the fact that the numbering
#		in the data file is 1 .. N, whereas here it is the
#		C convention, which is used 0 .. N-1.  Macros could
#		have been devised READ_FROM_R() and PRINT_TO_R().
#		As for now, the instances where such manipulations are
#		needed are flagged with 'C-numbering'.
#
#		Basically, as soon as possible everything is done
#		using the C convention.  Turning back to the R
#		convention is done only at the last minute at the time
#		of outputting something.
#
#		* Also, the matrices are flattened by R but in a
#              surprising way.  {{a, b}, {c, d}} becomes {a, c, b, d}.
#
#		* It is to be noted as well that the adjacency lists
#              describe incoming edges.  The neighbours of i
#              influence i!  And not: The neighbours of i are
#              influenced by i.  This is because the class of i is
#              influenced by its neighbours and not the opposite.
#
# ACKNOWLEDGMENTS: The author expresses his gratitude to Robert
#		Sochon (University of Sheffield, the UK) and William
#		Venables (CSIRO, Australia).

# FUNCTION: MMG.read.file(FILE-NAME)
#
#	This function reads the file FILE-NAME, which is supposed to
#	contain the data describing the adjacency list of the metabolic
#	network alongside the MS quantifications.  (NA should be used
#	for non-quantified proteins.)
#
#	The format of the file is:
#
#	    INDEX VALUE NEIGHBOUR1 W1 NEIGHBOUR2 W2 . . . NEIGHBOURn Wn
#
#	where:
#
#	    INDEX: is the index of the node
#	    VALUE: is a real value or NA (MS measurement)
#	    NEIGHBOURi: the indices of the neighbours (connected through a
#	                compound of weight Wi)
#
#	The value returned is a list of vectors representing the
#	adjacency matrix alongside the MS quantifications.  After
#
#	    dat <- MMG.read.file(FILE)
#
#	is a list with the following entries: adjacency.matrix, lengths, and
#	n.nodes.  Therefore:
#
#	    dat$adj[k, 2]: the MS quantification
#	    dat$adj[k, 2 i + 1]: the i-th neighbour
#	    dat$adj[k, 2 i + 2]: the i-th weight
#
#	    dat$len[k]: length of dat$adj[k]
#
#	    dat$n: number of nodes
#
#	Example:
#
#	    dat <- MMG.read.file("test.dat")

`MMG.read.file` <- function (file.name)
{
                                        # Read the lines
  conn <- file(file.name, "r");
  lins <- readLines(conn);
  close(conn);
  
                                        # Convert in numbers (list)
  f1 <- function (string) {
    as.numeric(strsplit(string, " +")[[1]])
  }
  adj.list <- lapply(as.list(lins), f1);
  n.nodes  <- length(adj.list);
  
                                        # Turn the list into a matrix
  f2 <- function (l) {
    length(l)
  }
  adj.len <- unlist(lapply(adj.list, f2));
  max.len <- max(adj.len);
  adj.mat <- matrix(0, nrow = n.nodes, ncol = max.len);

  for (i in 1:n.nodes) {
    for (j in 1:length(adj.list[[i]])) {
                                        # TODO: Invert the weights
					# TODO: Log(MS)
      adj.mat[i, j] <- adj.list[[i]][j]
    }
  }
  
  list(n.nodes = n.nodes, adjacency.matrix = adj.mat, lengths = adj.len)
}

# FUNCTION: MMG.compute(FILE-NAME, SIGMA, ALPHA, BURN-IN, STEPS)
#
#	Performs a Gibbs sampling based on the data obtained from
#	FILE-NAME and the model described in Sanguinetti & al,
#	Bioinfomatics 2008.  The following parameters can be set:
#
#	* SIGMA: Standard deviation of the unchanged distribution;
#
#	* ALPHA: Interplay between Standard Mixture Model and
#	Mixture Model on Graphs.  If ALPHA small, only the graph
#	structure counts, and if ALPHA -> +∞, one tends towards a
#	pure Mixture Model.
#
#	* BURN-IN and STEPS are the number of steps performed
#	by the Gibbs sampler.
#
# RETURNS: a list with the following fields,
#
#	* LUP
#	* LDOWN
#	* DATA

`MMG.compute` <- function(file.name, sigma = 0.3, alpha = 1,
                          burn.in = 1000, steps = 5000) {

# Read the data and create an array of samples
  dat     <- MMG.read.file(file.name);
  l       <- dim(dat$adj)[2];
  samples <- matrix(0, nrow = 3, ncol = dat$n.nodes);
  lup     <- rep(0.0, times = steps);
  ldown   <- rep(0.0, times = steps);

  r <- .C("MMG_compute",
  #                                                No. Meaning           (dims)    [C-equiv]
  # Data
    as.integer(c(dat$n.nodes, (l - 2) / 2)),     #  1. Dimensions N x L  (2)       [dims]
    as.integer((dat$len - 2) / 2),               #  2. Neighbour number  (N)       [nbnums]
    as.double  (dat$adj[, 2]),                   #  3. MS values         (N)       [msvals]
    as.integer (dat$adj[, seq(3, l, 2)]),        #  4. Neighbour indices (N x L)   [nbs]
    as.double  (dat$adj[, seq(4, l, 2)]),        #  5. Weights           (N x L)   [wts]
  # Parameters
    as.double(sigma),                            #  6.                             [sigma]
    as.double(alpha),                            #  7.                             [alpha]
    as.integer(burn.in),                         #  8.                             [burn_in]
    as.integer(steps),                           #  9.                             [steps]
  # Samples
    as.integer(samples),                         # 10.                             [samples]
    as.double(lup),                              # 11.                             [lup]
    as.double(ldown),                            # 12.                             [ldown]
  # Package information
    NAOK = TRUE,
    PACKAGE = "MMG"
  );

  p <- r[[10]] / steps;

# Print information
  cat("MMG Gibbs sampler (", steps, " steps)\n", sep = "");
  cat("Lambda down:     ", mean(r[[12]]), "; SD = ", sd(r[[12]]), "\n",
    sep = "");
  cat("Lambda up:       ", mean(r[[11]]), "; SD = ", sd(r[[11]]), "\n",
    sep = "");

# Shannon Entropy
  f <- function (x1, x2, x3) {
    ff(x1) + ff(x2) + ff(x3)
  }

  ff <- function (x) {
    if (x == 0.0)
      0
    else
      x * log(x)
  }

  e <- rep(0.0, times = dat$n.nodes);

  for (i in seq(1, dat$n.nodes)) {
    e[i] = -f(p[[(i - 1) * 3 + 1]],
              p[[(i - 1) * 3 + 2]],
              p[[(i - 1) * 3 + 3]])
  }

  cat("Shannon entropy: ", mean(e), "; SD = ", sd(e), "\n", sep = "");
  stem(e);

# Returns...
  res <- list(data = dat, samples = matrix(p, nrow = dat$n.nodes, ncol = 3, byrow = TRUE),
              lup = r[[11]], ldown = r[[12]],
              entropies = e);
  invisible(res);
}
