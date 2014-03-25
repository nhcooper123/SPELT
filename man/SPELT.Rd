\name{SPELT}
\alias{SPELT}
\title{SPELT}
\usage{
SPELT(phy, data, primary.variable, lag.variable, speciesnames, age.limit = NA,
  warn.dropped = TRUE, cut.off = 3)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}

  \item{data}{dataframe of class 'data.frame'}

  \item{node.list}{vector of node labels}

  \item{var1.col}{column number of variable 1}

  \item{var2.col}{column number of variable 2}

  \item{SPELT.data}{name of SPELT dataframe}

  \item{age.limit}{minimum age of branches to retain in
  SPELT. Default = NULL}

  \item{cut.off}{number of rows that must remain in the
  data for SPELT to continue. Default = 3}
}
\value{
An object of class 'SPELT' which contains the dataframe for
SPELT, the summary of the SPELT model, a list of species
dropped because they were not in both the phylogeny and
data, the primary variable, lag variable, and age limit
entered by the user, and the number of cherries.
}
\description{
First the function creates a dataframe via
\code{\link{get.SPELT.data}}. It then fits a model of
residuals from lag variables contrasts ~ primary variable
contrasts against divergence time for each species pair.
The optional argument \code{age.limit} allows the user to
remove species pairs with branches that are shorter (i.e.
younger) than a user defined age limit. The optional
argument \code{cut.off} allows users to decide how many
species pairs they need to trust the analyses (the default
is 3).
}
\examples{
data(shorebird)
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
summary(SPELT.results)
print(SPELT.results)
str(SPELT.results$data)

# With an age limit of 10 million years
SPELT.results10MY <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species", age.limit = 10)
summary(SPELT.results10MY)
print(SPELT.results10MY)
str(SPELT.results10MY$data)
}

