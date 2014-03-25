\name{build.SPELT.data}
\alias{build.SPELT.data}
\title{Build empty dataframe for SPELT functions}
\usage{
build.SPELT.data(phy)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}
}
\value{
An empty dataframe with the columns: "species1",
"species2", "species1.var1", "species2.var1",
"species1.var2", "species2.var2", "branch.length",
"contrast.var1", "contrast.var2", and "residuals".
}
\description{
Creates an empty dataframe to fill with data needed for
SPELT functions and adds column names. The rumber of rows
is equal to the number of cherries from
\code{\link{cherry.nodes}}.
}
\examples{
data(shorebird)
empty.dataset <- build.SPELT.data(shorebird.tree)
head(empty.dataset)
}
\seealso{
\code{\link{cherry.nodes}}
}

