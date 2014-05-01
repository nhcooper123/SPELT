\name{branch.length.pair}
\alias{branch.length.pair}
\title{Extract branch lengths leading to nodes of a phylogeny}
\usage{
branch.length.pair(phy, node.list)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}

  \item{node.list}{vector of node labels}
}
\value{
A vector of branch lengths
}
\description{
Finds the lengths of the branches leading to pair of
species subtending from a given set of nodes.
}
\examples{
data(shorebird, package="caper")
cherry.node.list <- cherry.nodes(shorebird.tree)
branch.length.pair(shorebird.tree, cherry.node.list)
}
\seealso{
\code{\link{cherry.nodes}}
}

