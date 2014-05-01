\name{fit.lag.model}
\alias{fit.lag.model}
\title{Fit lag residuals model}
\usage{
fit.lag.model(SPELT.data)
}
\arguments{
  \item{SPELT.data}{dataframe created by
  \code{\link{get.SPELT.data}}}
}
\value{
lm
}
\description{
Fits model of residuals from lag variables contrasts ~
primary varibale contrasts against divergence time for each
species pair.
}
\examples{
# see SPELT
}
\seealso{
\code{\link{get.SPELT.data}}
}

