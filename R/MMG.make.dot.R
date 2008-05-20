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

`MMG.make.dot` <- function (s, file.name, type = "DIRECTED",
selection, rem.loops = FALSE, weight.max = NA, weight.min = NA)
{
# C handles better numbers than strings
  tp <- if (type == "DIRECTED")
    1
  else {
    if (type == "UNDIRECTED")
      2
    else {
      warning("Invalid selection in MMG.make.dot(..., type = '%s', ...)\n", type);
      1
    }
  }

# Network topology
  dat     <- s$dat;
  l       <- dim(dat$adj)[2];
  samples <- s$samples;

  .C("MMG_make_dot",
          as.integer(c(dat$n.nodes, (l - 2) / 2)),     #  1. Dimensions N x L   (2)       [dims]
          as.integer(selection),                       #  2. Nodes to represent (K)       [select]
          as.integer(length(selection)),               #  3. Nodes to represent (K)       [selectn]
          as.character(file.name),                     #  4. File to output     (1)       [out]
          as.integer((dat$len - 2) / 2),               #  5. Neighbour number   (N)       [nbnums]
          as.integer(dat$adj[, seq(3, l, 2)]),         #  6. Neighbour indices  (N x L)   [nbs]
          as.logical(rem.loops),                       #  7. Remove loops       (1)       [rm_loops]
          as.integer(tp),                              #  8. DIRECTED or UNDIRECTED (1)   [type]
          as.double(samples),                          #  9. Samples            (N x 3)   [samples]
          PACKAGE = "MMG");

  return()
}
