\name{remove.young.branches}
\alias{remove.young.branches}
\title{Remove branches shorter than user defined age limit}
\usage{
remove.young.branches(SPELT.data, age.limit = NULL, cut.off = 3)
}
\arguments{
  \item{SPELT.data}{dataframe created by
  \code{\link{build.SPELT.data}} and filled with
  \code{\link{add.SPELT.data}}}

  \item{age.limit}{minimum age of branches to retain in
  SPELT. Default = NULL}

  \item{cut.off}{number of rows that must remain in the
  data for SPELT to continue. Default = 3}
}
\value{
A dataframe for SPELT where all nodes have branches > age
limit
}
\description{
Removes species pairs with branches that are shorter (i.e.
younger) than a user defined age limit. Practically this
removes rows from dataframe created by
\code{\link{build.SPELT.data}} and filled with
\code{\link{add.SPELT.data}}. Must be applied
\strong{before} adding the lag contrasts model with
\code{\link{get.SPELT.contrasts.data}}, as it will alter
the residuals. See \code{\link{get.SPELT.data}} for
examples.
}
\examples{
see get.SPELT.data
}

