\name{SPELT.summary.details}
\alias{SPELT.summary.details}
\title{Collate details for summary and plot outputs}
\usage{
SPELT.summary.details(SPELT.results)
}
\arguments{
  \item{SPELT.results}{object of class 'SPELT' from
  \code{\link{SPELT}}}
}
\value{
list of variable names and age limit
}
\description{
Collates names of variables and age limit for summary
functions.
}
\examples{
data(shorebird, package="caper")
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
summary(SPELT.results)
plot(SPELT.results)
}

