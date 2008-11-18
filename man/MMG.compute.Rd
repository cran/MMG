\name{MMG.compute}
\alias{MMG.compute}
\alias{MMG.cut.graph}
\alias{MMG.make.dot}
\title{Mixture Model on Graph (functions)}
\description{
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
  for visualisation using various softwares (including GraphViz
  \code{http://www.graphviz.org/} - see also the package
  \code{Rgraphviz}).
}
\usage{
r <- MMG.compute(file.name, data, sigma = 0.3, alpha = 1,
                 burn.in = 1000, steps = 5000)

s <- MMG.cut.graph(r, method = "THRESHOLD", select = "UP",
                   threshold = 0.2, descriptions = NA)

MMG.make.dot(s, file.name, type = "DIRECTED",
             selection, rem.loops = FALSE,
	     weight.max = NA, weight.min = NA)
}
\arguments{
  \item{r}{List returned by \code{MMG.compute}, cf. infra.}
  \item{s}{List returned by \code{MMG.cut.graph}, cf. infra.}

  \item{file.name}{\code{file.name} used in \code{MMG.compute} is the
             file containing the network alongside the (i.e.,
             proteomic) data.  See details for complementary details.
             Used in \code{MMG.make.dot}, it indicates the file that
             must be output.}

  \item{data}{Allows the user to override the data values using an
	additional file.  This may spare the user the inconvenience of
	re-generating network files all the time.}

  \item{sigma}{\code{sigma} is the standard deviation of the unchanged
               category (Gaussian).}

  \item{alpha}{\code{alpha} controls how important is the contribution
               of the network in the Mixture Model.  Zero means that
               the contribution is high whereas the model tends
               towards a normal mixture model when alpha becomes large.}

  \item{burn.in}{\code{burn.in} is the number of steps of burn-in of
                 the Gibbs sampler - if one feels that it can be useful.}

  \item{steps}{\code{steps} Number of steps that must be performed by
               the Biggs sampler.}

  \item{method}{\code{method} determines the method that will be used
                to identify the clusters of nodes that behave
                consistently, i.e., that belongs to a same class.
                Possible methods are "LIKELIEST", "THRESHOLD", and
                "ENTROPY".  See \code{threshold} and details.}

  \item{select}{\code{select} can be \code{"DOWN"},
                \code{"UNCHANGED"}, or \code{"UP"} depending on the
                class one desire to inspect.}

  \item{threshold}{\code{threshold} controls the stringency that is
                   used to determine whether a node belongs to the
                   class selected by \code{select}.  Depending on
                   \code{method}, the meaning of this parameter
                   changes.  See details.}

  \item{descriptions}{A file name.  The corresponding file must
                      contain the list of descriptions of the
                      network's nodes, one per line, in the same order
                      as \code{file.name}'s.  See details.}

  \item{type}{It can be \code{"DIRECTED"} or \code{"UNDIRECTED"}.
              Because the DOT fomat understands both types of graph.
              \code{"UNDIRECTED"} makes more sense if the original
              topology as specified by \code{file.name} in
              \code{MMG.compute} was undirected.  It can also make sense
              in the case of directed graphs if it makes them look better
              in GraphViz.  The graph is made undirected if necessary.}

  \item{selection}{This parameter is used to tell \code{MMG.make.dot}
                   which subset of the network must be represented.
                   It is generally computed using the value returned
                   by \code{MMG.cut.graph}.}

  \item{rem.loops}{A boolean value controlling whether loops
                   should be avoided.}

  \item{weight.max}{This parameter discards any edge that has a weight
                    exceeding it.  (Not implemented yet.)}

  \item{weight.min}{This parameter discards any edge that has a weight
                    lower than it.  (Not implemented yet.)}
}
\details{

  \code{MMG.compute} runs a Gibbs sampler and returns useful
    information regarding the posterior probabilities of belonging to
    such-and-such category (down-regulated, unchanged, and
    up-regulated), assuming the model described in Guido et al.,
    Bioinformatics (2008).

  The file containing the data must be made of lines having the following format:

  \preformatted{n value neighbour1 weight1 ... neighbourN weightN}

  \code{n} is the number of the line (this is required to make the
    file easier to read to a human being - it could be, in principle,
    any number); \code{value} is the logarithm of experimental
    measurement (0 meaning no change, and \code{NA} the value was not
    got), the base utilised does not matter very much but \code{sigma}
    should be chosen accordingly; the \code{neighbourI}s are the
    neighbours of the node \code{n} alongside the weights
    \code{weightI}s of the edges that connect it to them.  \code{N}
    needs not be the same all throughout the file.

  NB: It must be noted that contrarily to the weights used in Croes et
    al., JMB 2006, here a high weight implies high impact.  The weight
    used by Croes et al. must therefore be inverted before being used.

  \code{MMG.compute} prints some information regarding the parameters
  \itemize{
    \item{\code{lambda\_down} and \code{lambda\_up}}{Average and standard deviation}
    \item{Shannon entropy}{Average, standard deviation, and
          stem-and-leaves diagram.  This could give a quick flavour of
          how much uncertainty there is throughout the network.}
  }

  The file whose name is given by \code{descriptions} must contain
    lines, each of which, say the i-th line, describe the content of
    the node number i of the network.  For instance,

  \preformatted{aldehyde dehydrogenase
acetyl-coenzyme a synthetase
alcohol dehydrogenase
pyruvate dehydrogenase
dihydrolipoamide acetyltransferase
pyruvate kinase
...}

  The parameter \code{threshold} helps to select the nodes that belong
  to the desired category (parameter \code{select}).  Once the
  selection is carried out, \code{MMG.cut.graph} only has to generate
  the subgraph that contains the selection.  In the following T
  denotes the parameter \code{threshold}.

  \itemize{
    \item{\code{"LIKELIEST"}}{A node belongs to class C if and only if
      the probability pC of belonging to the class C is greater than
      that of belonging to the other classes D1 and D2 (pC>pD1 and
      pC>pD2), and if pC > T.}

    \item{\code{"THRESHOLD"}}{A node belongs to class C if and only if
      the probability pC of belonging to the class C is greater than
      that of belonging to the other classes D1 and D2 by at least T:
      pC>pD1+T and pC>pD2+T.}

  \item{\code{"ENTROPY"}}{A node belongs to class C if and only if
      the probability pC of belonging to the class C is the greatest
      (pC>pD1 and pC>pD2) and if the Shannon entropy is below T.}
  }
}
\value{
\code{MMG.compute} returns a list
\itemize{
  \item{data}{The data as read from the file \code{file.name}.  This
      is a list:
  \itemize{
    \item{n.nodes}{Number of nodes}
    \item{adjacency.matrix}{Contains a matrix each row of which is the list of values
      \preformatted{n value neighbour1 weight1 ... neighbourN weightN}
    }
    \item{lengths}{Number of fields of the rows}
  }}
  \item{samples}{The samples drawn from the Gibbs sampler (n x 3 matrix)}
  \item{lup}{The series of lambda\_up}
  \item{ldown}{The series of lambda\_down}
  \item{entropies}{List of the Shannon entropies}
}

\code{MMG.cut.graph} returns a list
\itemize{
  \item{components}{To which component belongs the nodes (0 means not selected)}
  \item{descriptions}{The vector of the descriptions as given in the
    file passed as \code{descriptions}}
}

\code{MMG.make.dot} does not return anything.
}
\references{
  Sanguinetti, Noirel, and Wright.,
  MMG: a probabilistic tool to identify submodules of metabolic
  pathways,
  Bioinformatics (2008)}
\author{Josselin Noirel, <j.noirel@sheffield.ac.uk>}
\note{
  Josselin Noirel,
  Departement of Chemical and Process Engineering,
  University of Sheffield, Mappin Street,
  Sheffield, S1 3JD - The United Kingdom,

  Based on an original implementation by Guido Sanguinetti
  <URL: http://sheffield.ac.uk/~guido/ >.
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ models }
