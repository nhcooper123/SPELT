\name{summary.SPELT}
\alias{summary.SPELT}
\title{Summary function for SPELT objects}
\usage{
\method{summary}{SPELT}(object, ...)
}
\arguments{
  \item{object}{object of class 'SPELT' from
    \code{\link{SPELT}}}
  \item{...}{Additional arguments to other methods (ignored)}  
}
\value{
\code{summary.lm} output
}
\description{
Summary output from a linear model of residuals from a
model of lag variable contrasts ~ primary variable
contrasts (through the origin), against divergence time for
all species pairs. Includes details of variables 
used and user defined \code{age.limit} for branch lengths
}
\examples{
data(shorebird, package = "caper")
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
summary(SPELT.results)
}

