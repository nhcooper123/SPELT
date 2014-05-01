\name{get.SPELT.data}
\alias{get.SPELT.data}
\title{Create SPELT dataset}
\usage{
get.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data,
  age.limit = NULL, cut.off)
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
A dataframe for SPELT with the columns: "species1",
"species2", "species1.var1", "species2.var1",
"species1.var2", "species2.var2", "branch.length",
"contrast.var1", "contrast.var2", and "residuals".
}
\description{
First the function creates an empty dataframe to fill with
data needed for SPELT functions and adds column names. The
rumber of rows is equal to the number of cherries. For each
cherry node, the function then fills the dataframe with
species names, primary and lag variable values, branch
lengths leading to the nodes, contrasts in primary and lag
variables, and residuals from a model of lag variable
contrasts ~ primary variable contrasts (forced through the
origin). Contrasts are calculated so that the primary
variable (variable 1) contrast is always positive, so the
function first finds which species of a pair has the
highest value for variable 1, then subtracts the other
species value from this. The order of subtraction is
retained when calculating the lag variable (variable 2)
contrast so these may be positive or negative. The optional
argument \code{age.limit} allows the user to removes
species pairs with branches that are shorter (i.e. younger)
than a user defined age limit. The optional argument
\code{cut.off} allows users to decide how many species
pairs they need to trust the analyses (the default is 3).
}
\examples{
data(shorebird, package="caper")
cherry.node.list <- cherry.nodes(shorebird.tree)
\dontrun{
# RGF: There is no way this can possily run - SPELT.data as argument
SPELT.data <- get.SPELT.data(shorebird.tree, shorebird.data,
                             cherry.node.list, 3, 5, SPELT.data)
str(SPELT.data)

# With an age limit of 10 million years
SPELT.data10MY <- get.SPELT.data(shorebird.tree, shorebird.data,
                             cherry.node.list, 3, 5, SPELT.data, age.limit = 10)
str(SPELT.data10MY)

# With an age limit of 15 million years
SPELT.data15MY <- get.SPELT.data(shorebird.tree, shorebird.data,
                             cherry.node.list, 3, 5, SPELT.data, age.limit = 15)
#This will not run because there are < 3 species pairs left in the data
}
}

