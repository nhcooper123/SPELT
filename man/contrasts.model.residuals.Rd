\name{contrasts.model.residuals}
\alias{contrasts.model.residuals}
\title{Fit lag contrasts model}
\usage{
contrasts.model.residuals(SPELT.data)
}
\arguments{
  \item{SPELT.data}{dataframe created by
  \code{\link{build.SPELT.data}} and filled with
  \code{\link{add.SPELT.data}} and
  \code{\link{get.raw.contrasts}}}
}
\value{
A vector of residuals from the model
}
\description{
Fits a linear model of lag variable contrasts ~ primary
variable contrasts, forced through the origin, and extracts
residuals from the model.
}
\examples{
see get.SPELT.data
}
\seealso{
\code{\link{get.SPELT.data}}
}

