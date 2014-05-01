\name{node.species1}
\alias{node.species1}
\title{Identify species coming from nodes of a phylogeny}
\usage{
node.species1(phy, node.list)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}

  \item{node.list}{vector of node labels}
}
\value{
A vector of species names
}
\description{
Identifies the species coming from a list of nodes then
selects either the first ones listed
(\code{\link{node.species1}}) or the second ones listed
(\code{\link{node.species2}}).
}
\examples{
data(shorebird, package="caper")
cherry.node.list <- cherry.nodes(shorebird.tree)
node.species1(shorebird.tree, cherry.node.list)
node.species2(shorebird.tree, cherry.node.list)
}
\seealso{
\code{\link{node.species2}}, \code{\link{cherry.nodes}}
}

