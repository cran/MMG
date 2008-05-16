\name{MMG-package}
\alias{MMG-package}
\alias{MMG}
\docType{package}
\title{Mixture Model on Graphs}
\description{
  MMG is a Mixture Model which can integrate the structure of a network
  in the statistical analysis of data.

  This implementation MMG assumes the existence of three classes of
  genes/proteins/etc. within a network: down-regulated, up-regulated,
  or unchanged.  The underlying data are log-ratios, hence relative,
  measurements.

  The package aims at identifying clusters of genes/proteins/etc. that
  behave consistently along the network's pathways.  This is done by
  implementing the Bayesian model described in the Reference
  Sanguinetti et al. (2008), by running a Gibbs sampler, and by
  cutting the graph.

  The Gibbs sampler is run by \code{MMG.compute}.  To cut the graph is
  done by \code{MMG.cut.graph}.  Finally, \code{MMG.make.dot}
  conveniently produces a DOT file, a file whose format can be used
  for visualation using various softwares (including GraphViz
  \code{http://www.graphviz.org/} - see also the package
  \code{Rgraphviz}).
}
\details{
  \tabular{ll}{
    Package: \tab MMG\cr
    Type:    \tab Package\cr
    Version: \tab 1.2.1\cr
    Date:    \tab 2008-05-16\cr
    License: \tab GPL-3+\cr
  }
}
\author{
  Josselin Noirel, based on an original implementation by Guido
  Sanguinetti (http://www.dcs.shef.ac.uk/~guido/).

  Maintainer: Josselin Noirel <j.noirel@sheffield.ac.uk>
}
\references{
  Sanguinetti, Noirel, and Wright.,
  MMG: a probabilistic tool to identify submodules of metabolic
  pathways,
  Bioinformatics (2008)}
\keyword{package}
\examples{
\dontrun{r <- MMG.compute(file.name = "NostocData/R_net.dat",
                 steps = 100000, burn.in = 1000,
                 sigma = 0.3, alpha = 1)}
\dontrun{n <- r$dat$n.nodes}
\dontrun{s <- MMG.cut.graph(r, descriptions = "NostocData/R_descr.dat",
                   method = "THRESHOLD", threshold = 0.15, select = "UP")}
\dontrun{l <- (1:n)[s$components != 0]}
\dontrun{MMG.make.dot(r, file.name = "nostoc.dot", selection = l, type = "UNDIRECTED",
             rem.loops = TRUE)}
}
