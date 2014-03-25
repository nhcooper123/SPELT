# Functions required for SPELT method for detecting evolutionary lag

#' Identify cherries
#'
#' Identifies all cherries, i.e. independent pairs of species from one node, 
#' in a phylogeny.
#'
#' @param phy phylogeny of class 'phylo'
#' @return A vector of node labels
#' @seealso \code{\link{total.nodes}}
#' @examples
#' data(shorebird)
#' cherry.nodes(shorebird.tree)
cherry.nodes <- function(phy) {
  names(which(table(phy$edge[, 1][phy$edge[, 2] <= total.nodes(phy)]) == 2))
}

#' Identify species coming from nodes of a phylogeny
#'
#' Identifies the species coming from a list of nodes then selects either the first
#' ones listed (\code{\link{node.species1}}) 
#' or the second ones listed (\code{\link{node.species2}}).
#'
#' @param phy phylogeny of class 'phylo'
#' @param node.list vector of node labels
#' @return A vector of species names
#' @seealso \code{\link{node.species2}}, \code{\link{cherry.nodes}}
#' @examples
#' data(shorebird)
#' cherry.node.list <- cherry.nodes(shorebird.tree)
#' node.species1(shorebird.tree, cherry.node.list)
#' node.species2(shorebird.tree, cherry.node.list)
node.species1 <- function(phy, node.list) {
  sapply(node.list, function(x) 
         phy$tip.label[phy$edge[, 2][which(phy$edge[, 1] == x)][1]])
}

#' Identify species coming from node
#'
#' Identifies the species coming from a given node then selects either the first
#' ones listed (\code{\link{node.species1}}) 
#' or the second ones listed (\code{\link{node.species2}}).
#'
#' @param phy phylogeny of class 'phylo'
#' @param node.list vector of node labels
#' @return A vector of species names
#' @seealso \code{\link{node.species1}}, \code{\link{cherry.nodes}}
#' @examples
#' data(shorebird)
#' cherry.node.list <- cherry.nodes(shorebird.tree)
#' node.species1(shorebird.tree, cherry.node.list)
#' node.species2(shorebird.tree, cherry.node.list)
node.species2 <- function(phy, node.list) {
  sapply(node.list, function(x) 
         phy$tip.label[phy$edge[, 2][which(phy$edge[, 1] == x)][2]])
}

#' Extract branch lengths leading to nodes of a phylogeny
#'
#' Finds the lengths of the branches leading to pair of species subtending from 
#' a given set of nodes. 
#'
#' @param phy phylogeny of class 'phylo'
#' @param node.list vector of node labels
#' @return A vector of branch lengths
#' @seealso \code{\link{cherry.nodes}}
#' @examples
#' data(shorebird)
#' cherry.node.list <- cherry.nodes(shorebird.tree)
#' branch.length.pair(shorebird.tree, cherry.node.list)
branch.length.pair <- function(phy, node.list) {
  sapply(node.list, function(x) 
         phy$edge.length[which(phy$edge[,1] == x)][1])
}

#' Build empty dataframe for SPELT functions
#'
#' Creates an empty dataframe to fill with data needed for SPELT functions
#' and adds column names. The rumber of rows is equal to the number of cherries from \code{\link{cherry.nodes}}.
#'
#' @param phy phylogeny of class 'phylo'
#' @return An empty dataframe with the columns: "species1", "species2", "species1.var1", 
#' "species2.var1", "species1.var2", "species2.var2", "branch.length", 
#' "contrast.var1", "contrast.var2", and "residuals".
#' @seealso \code{\link{cherry.nodes}}
#' @examples
#' data(shorebird)
#' empty.dataset <- build.SPELT.data(shorebird.tree)
#' head(empty.dataset)
build.SPELT.data <- function(phy) {
  SPELT.data <- data.frame(array(dim = c(length(cherry.nodes(phy)), 10)))
  colnames(SPELT.data) <- c("species1", "species2", "species1.var1", 
                            "species2.var1", "species1.var2",
                            "species2.var2", "branch.length", 
                            "contrast.var1", "contrast.var2", 
                            "residuals")
  return(SPELT.data)
}

#' Add species and variable data into SPELT dataframe
#'
#' Fills empty dataframe created by  \code{\link{build.SPELT.data}} 
#' with species names, variable values and branch lengths. See \code{\link{get.SPELT.data}} for examples.
#'
#' @param phy phylogeny of class 'phylo'
#' @param data dataset of class 'data.frame'
#' @param node.list vector of node labels
#' @param var1.col column number of variable 1
#' @param var2.col column number of variable 2
#' @param SPELT.data empty dataframe created by \code{\link{build.SPELT.data}}
#' @return A mostly completed dataframe for SPELT
#' @seealso \code{\link{build.SPELT.data}}, \code{\link{get.raw.contrasts}}, 
#' \code{\link{node.species1}}, \code{\link{node.species2}}, \code{\link{branch.length.pair}},
#' \code{\link{get.SPELT.data}}
#' @examples
#' see get.SPELT.data
add.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, SPELT.data) {
  SPELT.data$species1 <- node.species1(phy,node.list)
  SPELT.data$species2 <- node.species2(phy,node.list)
  SPELT.data$species1.var1 <- get.data(data, var1.col, SPELT.data$species1)
  SPELT.data$species2.var1 <- get.data(data, var1.col, SPELT.data$species2)
  SPELT.data$species1.var2 <- get.data(data, var2.col, SPELT.data$species1)
  SPELT.data$species2.var2 <- get.data(data, var2.col, SPELT.data$species2)
  SPELT.data$branch.length <- branch.length.pair(phy, node.list)
  return(SPELT.data)
}

#' Calculate contrasts in primary and lag variables
#'
#' Adds contrasts to dataset created by \code{\link{build.SPELT.data}} and filled with
#' \code{\link{add.SPELT.data}}. Contrasts are first calculated using data already
#' entered into the dataset. The primary variable (variable 1) contrast must always be
#' positive, so the function first finds which species of a pair has the highest value
#' for variable 1, then subtracts the other species value from this. The order of 
#' subtraction is retained when calculating the lag variable (variable 2) contrast
#' so these may be positive or negative. See \code{\link{get.SPELT.data}} for examples.
#'
#' @param SPELT.data dataframe created by \code{\link{build.SPELT.data}} and filled with \code{\link{add.SPELT.data}} 
#' @return A mostly completed dataframe for SPELT
#' @seealso \code{\link{build.SPELT.data}}, \code{\link{add.SPELT.data}}, \code{\link{get.SPELT.data}}
#' @examples
#' see get.SPELT.data
get.raw.contrasts <- function(SPELT.data) { 
  for(i in seq_along(SPELT.data$species1)) {
    if (SPELT.data$species1.var1[i] >= SPELT.data$species2.var1[i]) {
      SPELT.data$contrast.var1[i] <- SPELT.data$species1.var1[i] - SPELT.data$species2.var1[i]
      SPELT.data$contrast.var2[i] <- SPELT.data$species1.var2[i] - SPELT.data$species2.var2[i]
    } else {
      SPELT.data$contrast.var1[i] <- SPELT.data$species2.var1[i] - SPELT.data$species1.var1[i]
      SPELT.data$contrast.var2[i] <- SPELT.data$species2.var2[i] - SPELT.data$species1.var2[i]
    }
  }  
  return(SPELT.data)
}

#' Fit lag contrasts model
#'
#' Fits a linear model of lag variable contrasts ~ primary variable contrasts, 
#' forced through the origin, and extracts residuals from the model. 
#'
#' @param SPELT.data dataframe created by \code{\link{build.SPELT.data}} 
#' and filled with \code{\link{add.SPELT.data}} and \code{\link{get.raw.contrasts}}
#' @return A vector of residuals from the model
#' @seealso \code{\link{get.SPELT.data}}
#' @examples
#' see get.SPELT.data
contrasts.model.residuals <- function(SPELT.data) {
  model <- lm(SPELT.data$contrast.var2 ~ SPELT.data$contrast.var1 - 1)
  residuals <- unclass(model)$residuals
  return(residuals)
}

#' Add contrasts residuals to SPELT dataset
#'
#' Adds residuals from lag contrasts model to dataset created by \code{\link{build.SPELT.data}} 
#' and filled with \code{\link{add.SPELT.data}} and \code{\link{get.raw.contrasts}}. 
#' See \code{\link{get.SPELT.data}} for examples.
#'
#' @param SPELT.data dataset created by \code{\link{build.SPELT.data}} 
#' and filled with \code{\link{add.SPELT.data}} and \code{\link{get.raw.contrasts}}
#' @return A completed dataframe for SPELT
#' @seealso \code{\link{build.SPELT.data}}, \code{\link{add.SPELT.data}}, \code{\link{get.raw.contrasts}},
#' \code{\link{get.SPELT.data}}
#' @examples
#' see get.SPELT.data
add.SPELT.contrasts.data <- function(SPELT.data) {
  SPELT.data$residuals <- contrasts.model.residuals(SPELT.data)
  return(SPELT.data)
}

#' Remove branches shorter than user defined age limit
#'
#' Removes species pairs with branches that are shorter (i.e. younger) than a user defined age limit.
#' Practically this removes rows from dataframe created by \code{\link{build.SPELT.data}} 
#' and filled with \code{\link{add.SPELT.data}}. Must be applied \strong{before} adding the lag
#' contrasts model with \code{\link{get.SPELT.contrasts.data}}, as it will alter the residuals.
#' See \code{\link{get.SPELT.data}} for examples.
#'
#' @param SPELT.data dataframe created by \code{\link{build.SPELT.data}} 
#' and filled with \code{\link{add.SPELT.data}}
#' @param age.limit minimum age of branches to retain in SPELT. Default = NULL
#' @param cut.off number of rows that must remain in the data for SPELT to continue. Default = 3
#' @return A dataframe for SPELT where all nodes have branches > age limit
#' @examples
#' see get.SPELT.data
remove.young.branches <- function(SPELT.data, age.limit = NULL, cut.off = 3) {
  if (!is.null(age.limit)) {
    SPELT.data <- SPELT.data[-(c(which(SPELT.data$branch.length < age.limit))), ]
  }
    if (nrow(SPELT.data) < cut.off) {
      stop(paste("< ", cut.off, " branches longer than age limit"))
    }
  return(SPELT.data) 
}

#' Create SPELT dataset
#'
#' First the function creates an empty dataframe to fill with data needed for SPELT functions
#' and adds column names. The rumber of rows is equal to the number of cherries. For each cherry node, the function then
#' fills the dataframe with species names, primary and lag variable values, branch lengths leading to the nodes, 
#' contrasts in primary and lag variables, and
#' residuals from a model of lag variable contrasts ~ primary variable contrasts (forced through the origin).
#' Contrasts are calculated so that the primary variable (variable 1) contrast is always
#' positive, so the function first finds which species of a pair has the highest value
#' for variable 1, then subtracts the other species value from this. The order of 
#' subtraction is retained when calculating the lag variable (variable 2) contrast
#' so these may be positive or negative. The optional argument \code{age.limit} allows the user to removes species pairs 
#' with branches that are shorter (i.e. younger) than a user defined age limit. The optional argument \code{cut.off}
#' allows users to decide how many species pairs they need to trust the analyses (the default is 3).
#'
#' @param phy phylogeny of class 'phylo'
#' @param data dataframe of class 'data.frame'
#' @param node.list vector of node labels
#' @param var1.col column number of variable 1
#' @param var2.col column number of variable 2
#' @param SPELT.data name of SPELT dataframe
#' @param age.limit minimum age of branches to retain in SPELT. Default = NULL
#' @param cut.off number of rows that must remain in the data for SPELT to continue. Default = 3
#' @return A dataframe for SPELT with the columns: "species1", "species2", "species1.var1", 
#' "species2.var1", "species1.var2", "species2.var2", "branch.length", 
#' "contrast.var1", "contrast.var2", and "residuals".
#'  @examples
#' data(shorebird)
#' cherry.node.list <- cherry.nodes(shorebird.tree)
#' SPELT.data <- get.SPELT.data(shorebird.tree, shorebird.data, 
#'                              cherry.node.list, 3, 5, SPELT.data)
#' str(SPELT.data)
#' 
#' # With an age limit of 10 million years
#' SPELT.data10MY <- get.SPELT.data(shorebird.tree, shorebird.data, 
#'                              cherry.node.list, 3, 5, SPELT.data, age.limit = 10)
#' str(SPELT.data10MY)
#'
#' # With an age limit of 15 million years
#' SPELT.data15MY <- get.SPELT.data(shorebird.tree, shorebird.data, 
#'                              cherry.node.list, 3, 5, SPELT.data, age.limit = 15)
#' #This will not run because there are < 3 species pairs left in the data
get.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, 
                           SPELT.data, age.limit = NULL, cut.off) {
  SPELT.data <- build.SPELT.data(phy)
  SPELT.data <- add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)
  SPELT.data <- get.raw.contrasts(SPELT.data)
  SPELT.data <- remove.young.branches(SPELT.data, age.limit, cut.off)
  SPELT.data <- add.SPELT.contrasts.data(SPELT.data)
  return(SPELT.data)
}

#' Fit lag residuals model
#'
#' Fits model of residuals from lag variables contrasts ~ primary varibale contrasts 
#' against divergence time for each species pair.
#'
#' @param SPELT.data dataframe created by \code{\link{get.SPELT.data}}
#' @return lm
#' @seealso \code{\link{get.SPELT.data}}
#' @examples
#' see SPELT
fit.lag.model <- function(SPELT.data) {
  lag.model <- lm(SPELT.data$residuals ~ SPELT.data$branch.length)
}

#' SPELT function
#'
#' First the function creates a dataframe via \code{\link{get.SPELT.data}}. It then fits a model of 
#' residuals from lag variables contrasts ~ primary variable contrasts 
#' against divergence time for each species pair.
#' The optional argument \code{age.limit} allows the user to removes species pairs 
#' with branches that are shorter (i.e. younger) than a user defined age limit. 
#' The optional argument \code{cut.off}
#' allows users to decide how many species pairs they need to trust the analyses (the default is 3).
#'
#' @param phy phylogeny of class 'phylo'
#' @param data dataframe of class 'data.frame'
#' @param node.list vector of node labels
#' @param var1.col column number of variable 1
#' @param var2.col column number of variable 2
#' @param SPELT.data name of SPELT dataframe
#' @param age.limit minimum age of branches to retain in SPELT. Default = NULL
#' @param cut.off number of rows that must remain in the data for SPELT to continue. Default = 3
#' @return An object of class 'SPELT' which contains the dataframe for SPELT, the summary of the
#' SPELT model, a list of species dropped becasue they were not in both the phylogeny and data, 
#' the primary variable, lag variable, and age limit entered by the user, and the number of cherries.
#' @examples
#' data(shorebird)
#' SPELT.results <- SPELT(shorebird.tree, shorebird.data, 
#'                          "F.Mass", "Egg.Mass", "Species")
#' summary(SPELT.results)
#' print(SPELT.results)
#' str(SPELT.results$data)
#'
#' # With an age limit of 10 million years
#' SPELT.results10MY <- SPELT(shorebird.tree, shorebird.data, 
#'                          "F.Mass", "Egg.Mass", "Species", age.limit = 10)
#' summary(SPELT.results10MY)
#' print(SPELT.results10MY)
#' str(SPELT.results10MY$data)
SPELT <- function(phy, data, primary.variable, lag.variable, speciesnames, 
                  age.limit = NA, warn.dropped = TRUE, cut.off = 3) {

  if (!is.data.frame(data)) 
    stop("'data' must be an object of class 'data.frame'")

  if (!inherits(phy, "phylo")) 
        stop("'phy' must be an object of class 'phylo'")

  # Define variables
  var1.col <- column.ID(data, primary.variable)
  if (length(var1.col) == 0)
    stop("Primary variable not found in data")  
  
  var2.col <- column.ID(data, lag.variable)
  if (length(var2.col) == 0)
    stop("Lag variable not found in data")
  
  speciesnames.col <- column.ID(data, speciesnames)
  if (length(speciesnames.col) == 0)
    stop("Species names not found in data")

  # Ensure phylogeny has no polytomies
  phy <- multi2di(phy)
  
  # Tidy up data and tree
  data <- remove.incomplete.data(data, var1.col, var2.col)
  data <- remove.missing.species.data(phy, data, speciesnames.col)
  if (nrow(data) < cut.off)
    stop(paste("< ", cut.off, " species have data for both variables and are in the phylogeny"))
  phy <- remove.missing.species.tree(phy, data, speciesnames.col)

  if (warn.dropped) {
    tree.not.data <- id.missing.tree(phy, data, speciesnames.col)
    data.not.tree <- id.missing.data(phy, data, speciesnames.col)
  }

  # Identify cherries
  node.list <- cherry.nodes(phy)
  
  # Collate data required for SPELT analyses
  SPELT.data <- get.SPELT.data(phy, data, node.list, 
                               var1.col, var2.col, SPELT.data, age.limit, cut.off)

  # Fit SPELT model 
  SPELT.model <- fit.lag.model(SPELT.data)
  
  # Outputs
  SPELT.results <- list(summary = summary(SPELT.model), variables = list(primary.variable = primary.variable, 
                        lag.variable = lag.variable), data = SPELT.data, age.limit = age.limit,
                        dropped = list(tree.not.data = tree.not.data, data.not.tree = data.not.tree), 
                        Nnodes = nrow(SPELT.data)) 

  class(SPELT.results) <- "SPELT"
  return(SPELT.results)
}

#' Collate details for summary and plot outputs
#'
#' Collates names of variables and age limit for summary functions.
#'
#' @param SPELT.results object of class 'SPELT' from \code{\link{SPELT}}
#' @return list of variable names and age limit
#' @examples
#' data(shorebird)
#' SPELT.results <- SPELT(shorebird.tree, shorebird.data, 
#'                          "F.Mass", "Egg.Mass", "Species")
#' summary(SPELT.results)
#' plot(SPELT.results)
SPELT.summary.details <- function(SPELT.results) {
  details <- paste("primary variable = ", SPELT.results$variables$primary.variable,
                      ", lag variable = ", SPELT.results$variables$lag.variable, sep = "")
  if (!is.null(SPELT.results$age.limit)) {
    age.limit <- paste("age limit = ", SPELT.results$age.limit, sep = "")
  } else {
    age.limit <- "age limit = NULL"
  }
  return(list(details, age.limit))
}

#' Plotting function for SPELT objects
#'
#' Plots residuals from a model of lag variable contrasts ~ primary variable contrasts (through the origin),
#' against divergence time for all species pairs.
#'
#' @param SPELT.results object of class 'SPELT' from \code{\link{SPELT}}
#' @return plot with model as a solid line, and the 0,0 line as a dashed line
#' with details of variable and age.limit
#' @examples
#' data(shorebird)
#' SPELT.results <- SPELT(shorebird.tree, shorebird.data, 
#'                          "F.Mass", "Egg.Mass", "Species")
#' plot(SPELT.results)
plot.SPELT <- function(SPELT.results) {
  par(bty = "l")
  plot(SPELT.results$data$residuals ~ SPELT.results$data$branch.length, 
       xlab = paste("divergence time (", SPELT.summary.details(SPELT.results)[[2]],")", sep = ""),
       ylab = "residuals", cex.main = 0.8, pch = 16, las = 1,
       main = paste("SPELT results: ", SPELT.summary.details(SPELT.results)[[1]], sep = ""))
  abline(fit.lag.model(SPELT.results$data))
  abline(0,0,lty = 2)
}  

#' Summary function for SPELT objects
#'
#' Summary output from a linear model of residuals from a model of lag variable contrasts ~ primary variable contrasts (through the origin),
#' against divergence time for all species pairs.
#'
#' @param SPELT.results object of class 'SPELT' from \code{\link{SPELT}}
#' @return summary.lm output with details of variable and age.limit
#' @examples
#' data(shorebird)
#' SPELT.results <- SPELT(shorebird.tree, shorebird.data, 
#'                          "F.Mass", "Egg.Mass", "Species")
#' plot(SPELT.results)
summary.SPELT <- function(SPELT.results) {
  cat("\nSPELT Details:\n", SPELT.summary.details(SPELT.results)[[1]])
  cat("\n",SPELT.summary.details(SPELT.results)[[2]], "\n")
  print(SPELT.results$summary)
}