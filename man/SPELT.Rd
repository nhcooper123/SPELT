\name{SPELT}
\alias{SPELT}
\title{SPELT function}
\usage{
SPELT(phy, data, primary.variable, lag.variable, speciesnames, age.limit = NULL,
  warn.dropped = TRUE, cut.off = 3)
}
\arguments{
  \item{phy}{phylogeny of class 'phylo'}

  \item{data}{dataframe of class 'data.frame'}

  \item{primary.variable}{column name of the variable thought to drive
  evolution of the lag variable}

  \item{lag.variable}{column name of the variable thought to lag
  behind the primary variable}

  \item{speciesnames}{column name of the column containing species names}

  \item{age.limit}{minimum age of branches to retain in
  SPELT. Default = NULL}

  \item{warn.dropped}{If TRUE the function prints a list of species
  found in the phylogeny but not the data and vice versa. Default = TRUE}

  \item{cut.off}{number of rows that must remain in the
  data for SPELT to continue. Default = 3}
}
\value{
An object of class 'SPELT' which contains the dataframe for
SPELT, the summary of the SPELT model, a list of species
dropped becasue they were not in both the phylogeny and
data, the primary variable, lag variable, and age limit
entered by the user, and the number of cherries.
}
\description{
The function fits a model of
residuals from lag variable contrasts ~ primary variable
contrasts against divergence time for each species pair.
The optional argument \code{age.limit} allows the user to
removes species pairs with branches that are shorter (i.e.
younger) than a user defined age limit. The optional
argument \code{cut.off} allows users to decide how many
species pairs they need to trust the analyses (the default
is 3).
}
\examples{
data(shorebird, package = "caper")
SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
summary(SPELT.results)
plot(SPELT.results)
str(SPELT.results$data)

# With an age limit of 10 million years
SPELT.results10MY <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species", age.limit = 10)
summary(SPELT.results10MY)
plot(SPELT.results10MY)
str(SPELT.results10MY$data)
}

