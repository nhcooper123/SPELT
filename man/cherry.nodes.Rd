\name{cherry.nodes}
\alias{cherry.nodes}
\title{Identify cherries}
\usage{
cherry.nodes(phy)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}
}
\value{
A vector of node labels
}
\description{
Identifies all cherries, i.e. independent pairs of species
from one node, in a phylogeny.
}
\examples{
data(shorebird)
cherry.nodes(shorebird.tree)
}
\seealso{
\code{\link{total.nodes}}
}

