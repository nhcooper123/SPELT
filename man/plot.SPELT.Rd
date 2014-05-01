\name{plot.SPELT}
\alias{plot.SPELT}
\title{Plotting function for SPELT objects}
\usage{
\method{plot}{SPELT}(x, ...)
}
\arguments{
  \item{x}{object of class 'SPELT' from
    \code{\link{SPELT}}}
  \item{...}{Additional arguments to other methods (ignored)}
}
\value{
plot with model as a solid line, and the 0,0 line as a
dashed line with details of variable and age.limit
}
\description{
Plots residuals from a model of lag variable contrasts ~
primary variable contrasts (through the origin), against
divergence time for all species pairs.
}
\examples{
data(shorebird, package="caper")
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
plot(SPELT.results)
}

