\name{add.SPELT.contrasts.data}
\alias{add.SPELT.contrasts.data}
\title{Add contrasts residuals to SPELT dataset}
\usage{
add.SPELT.contrasts.data(SPELT.data)
}
\arguments{
  \item{SPELT.data}{dataset created by
  \code{\link{build.SPELT.data}} and filled with
  \code{\link{add.SPELT.data}} and
  \code{\link{get.raw.contrasts}}}
}
\value{
A completed dataframe for SPELT
}
\description{
Adds residuals from lag contrasts model to dataset created
by \code{\link{build.SPELT.data}} and filled with
\code{\link{add.SPELT.data}} and
\code{\link{get.raw.contrasts}}. See
\code{\link{get.SPELT.data}} for examples.
}
\examples{
see get.SPELT.data
}
\seealso{
\code{\link{build.SPELT.data}},
\code{\link{add.SPELT.data}},
\code{\link{get.raw.contrasts}},
\code{\link{get.SPELT.data}}
}

