\name{get.raw.contrasts}
\alias{get.raw.contrasts}
\title{Calculate contrasts in primary and lag variables}
\usage{
get.raw.contrasts(SPELT.data)
}
\arguments{
  \item{SPELT.data}{dataframe created by
  \code{\link{build.SPELT.data}} and filled with
  \code{\link{add.SPELT.data}}}
}
\value{
A mostly completed dataframe for SPELT
}
\description{
Adds contrasts to dataset created by
\code{\link{build.SPELT.data}} and filled with
\code{\link{add.SPELT.data}}. Contrasts are first
calculated using data already entered into the dataset. The
primary variable (variable 1) contrast must always be
positive, so the function first finds which species of a
pair has the highest value for variable 1, then subtracts
the other species value from this. The order of subtraction
is retained when calculating the lag variable (variable 2)
contrast so these may be positive or negative. See
\code{\link{get.SPELT.data}} for examples.
}
\examples{
see get.SPELT.data
}
\seealso{
\code{\link{build.SPELT.data}},
\code{\link{add.SPELT.data}}, \code{\link{get.SPELT.data}}
}

