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
summary.lm output with details of variable and age.limit
}
\description{
Summary output from a linear model of residuals from a
model of lag variable contrasts ~ primary variable
contrasts (through the origin), against divergence time for
all species pairs.
}
\examples{
data(shorebird, package="caper")
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
plot(SPELT.results)
}

