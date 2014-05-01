\name{add.SPELT.data}
\alias{add.SPELT.data}
\title{Add species and variable data into SPELT dataframe}
\usage{
add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}

  \item{data}{dataset of class 'data.frame'}

  \item{node.list}{vector of node labels}

  \item{var1.col}{column number of variable 1}

  \item{var2.col}{column number of variable 2}

  \item{SPELT.data}{empty dataframe created by
  \code{\link{build.SPELT.data}}}
}
\value{
A mostly completed dataframe for SPELT
}
\description{
Fills empty dataframe created by
\code{\link{build.SPELT.data}} with species names, variable
values and branch lengths. See \code{\link{get.SPELT.data}}
for examples.
}
\examples{
# see get.SPELT.data
}
\seealso{
\code{\link{build.SPELT.data}},
\code{\link{get.raw.contrasts}},
\code{\link{node.species1}}, \code{\link{node.species2}},
\code{\link{branch.length.pair}},
\code{\link{get.SPELT.data}}
}

