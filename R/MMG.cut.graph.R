# MODULE:      MMG
#
# FILE:	MMG_make_dot
#
# DESCRIPTION:	Outputs a DOT file with the network in it.
#
# VERSION:     1.2.2
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

# PARAMETERS
#
#	* data, as returned by MMG.compute()
#
#	* method = "LIKELIEST" | "THRESHOLD" | "ENTROPY"
#
#	* select = "UP" | "DOWN" | "UNCHANGED"

# METHODS
#
# method = "LIKELIEST", threshold = THRESHOLD, select = CLASS
#
#	We select the nodes such that P(CLASS) > THRESHOLD
#
# method = "THRESHOLD", threshold = THRESHOLD, select = CLASS
#
#	We select the nodes such that P(CLASS) > P(other classes) + THRESHOLD
#
# method = "ENTROPY", threshold = THRESHOLD, select = CLASS
#
#	We select the nodes such that P(CLASS) > P(other classes) && Shannon_Entropy > THRESHOLD

`MMG.cut.graph` <- function (data, method = "THRESHOLD", select = "UP", threshold = 0.2, descriptions = NA)
{
# C handles better numbers than strings
  sel <- if (select == "DOWN")
    1
  else {
    if (select == "UNCHANGED")
      2
    else {
      if (select == "UP")
        3
      else {
        warning("Invalid selection in MMG.cut.graph(..., select = '%s', ...)\n", select);
        3
      }
    }
  }

# C handles better numbers than strings
  met <- if (method == "LIKELIEST")
    1
  else {
    if (method == "THRESHOLD")
      2
    else {
      if (method == "ENTROPY")
        3
      else {
        warning("Invalid method in MMG.cut.graph(..., method = '%s', ...)\n", method);
        3
      }
    }
  }

# Samples from the Gibbs sampler
  samples <- data$samples;

# Network topology
  dat     <- data$dat;
  l       <- dim(dat$adj)[2];
  categ   <- rep(0, times = dat$n.nodes); # Which connected component is belongs to?
  ncateg  <- 0;

# Read the descriptions
  if (! is.na(descriptions)) {
    conn <- file(descriptions);
    dsc  <- readLines(conn);
    close(conn);

    if (length(dsc) != dat$n.nodes)
      warning("The length of the data and of the descriptions do not match!\n");
  }
  else
    dsc  <- NA;

  r <- .C("MMG_cut_graph",
          as.integer(c(dat$n.nodes, (l - 2) / 2)),     #  1. Dimensions N x L  (2)       [dims]
          as.integer((dat$len - 2) / 2),               #  2. Neighbour number  (N)       [nbnums]
          as.double  (dat$adj[, 2]),                   #  3. MS values         (N)       [msvals]
          as.integer (dat$adj[, seq(3, l, 2)]),        #  4. Neighbour indices (N x L)   [nbs]
          as.double  (dat$adj[, seq(4, l, 2)]),        #  5. Weights           (N x L)   [wts]
          as.double  (samples),                        #  6. Samples           (N x 3)   [samples]
          as.integer(met),                             #  7. METHOD     1, 2, or 3       [method]
          as.integer(sel),                             #  8. SELECTION  1, 2, or 3       [select]
          as.double(threshold),                        #  9. THRESHOLD                   [thres]
          as.integer(categ),                           # 10. categories        (N)       [categ]
          as.integer(ncateg),                          # 11. categories num    (1)       [ncateg]
          NAOK = TRUE,
          PACKAGE = "MMG");

  if (r[[11]] > 1) {
    for (cat in 1:(r[[11]] - 1)) {
      cat("COMPONENT ", cat, "\n", sep = "");
      for (i in 1:dat$n.nodes) {
        if (r[[10]][i] == cat) {
          ss <- if (is.na(dsc[1]))
            sprintf("  * Node %3d %.3lf [%1.3lf %1.3lf %1.3lf]\n",    i, dat$adjacency[i, 2], samples[i, 1], samples[i, 2], samples[i,3])
          else
            sprintf("  * Node %3d %.3lf [%1.3lf %1.3lf %1.3lf] %s\n", i, dat$adjacency[i, 2], samples[i, 1], samples[i, 2], samples[i,3], dsc[i]);
          cat(ss);
        }
      }
      cat("\n");
    }
  }

  s <- list(components = r[[10]],
            descriptions = dsc);
  invisible(s)
}
