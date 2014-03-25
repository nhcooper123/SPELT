# Functions required for SPELT method for detecting evolutionary lag

# Identify cherries (independent pairs of species from one node)
cherry.nodes <- function(phy) {
  names(which(table(phy$edge[, 1][phy$edge[, 2] <= total.nodes(phy)]) == 2))
}

# Identify 1st species coming from node
node.species1 <- function(phy, node.list) {
  sapply(node.list, function(x) 
         phy$tip.label[phy$edge[, 2][which(phy$edge[, 1] == x)][1]])
}

# Identify 2nd species coming from node
node.species2 <- function(phy, node.list) {
  sapply(node.list, function(x) 
         phy$tip.label[phy$edge[, 2][which(phy$edge[, 1] == x)][2]])
}

# Extract branch length for contrast between two species
branch.length.pair <- function(node.list) {
  sapply(node.list, function(x) 
         phy$edge.length[which(phy$edge[,1] == x)][1])
}

# Build empty dataset for SPELT
build.SPELT.data <- function(phy) {
  SPELT.data <- data.frame(array(dim = c(length(cherry.nodes(phy)), 10)))
  colnames(SPELT.data) <- c("species1", "species2", "species1.var1", 
                            "species2.var1", "species1.var2",
                            "species2.var2", "branch.length", 
                            "contrast.var1", "contrast.var2", 
                            "residuals")
  return(SPELT.data)
}

# Add species and variable data into SPELT dataset
add.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, SPELT.data) {
  SPELT.data$species1 <- node.species1(phy,node.list)
  SPELT.data$species2 <- node.species2(phy,node.list)
  SPELT.data$species1.var1 <- get.data(data, var1.col, SPELT.data$species1)
  SPELT.data$species2.var1 <- get.data(data, var1.col, SPELT.data$species2)
  SPELT.data$species1.var2 <- get.data(data, var2.col, SPELT.data$species1)
  SPELT.data$species2.var2 <- get.data(data, var2.col, SPELT.data$species2)
  SPELT.data$branch.length <- branch.length.pair(node.list)
  return(SPELT.data)
}

# Calculate contrasts (primary variable contrast is always positive)
# And add into SPELT dataset
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

# Fit model of lag variable contrasts ~ primary variable contrasts (forced through the origin)
contrasts.model.residuals <- function(SPELT.data) {
  model <- lm(SPELT.data$contrast.var2 ~ SPELT.data$contrast.var1 - 1)
  residuals <- unclass(model)$residuals
  return(residuals)
}

# Add contrasts residuals to SPELT dataset
add.SPELT.contrasts.data <- function(SPELT.data) {
  SPELT.data$residuals <- contrasts.model.residuals(SPELT.data)
  return(SPELT.data)
}

# Remove branches shorter than user defined age limit
remove.young.branches <- function(SPELT.data, age.limit = NULL, cut.off = 3) {
  if (!is.null(age.limit)) {
    SPELT.data <- SPELT.data[-(c(which(SPELT.data$branch.length < age.limit))), ]
  }
    if (nrow(SPELT.data) < cut.off) {
      stop(paste("< ", cut.off, " branches longer than age limit"))
    }
  return(SPELT.data) 
}

# Combining all functions for data collation
get.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, 
                           SPELT.data, age.limit = NULL, cut.off) {
  SPELT.data <- build.SPELT.data(phy)
  SPELT.data <- add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)
  SPELT.data <- get.raw.contrasts(SPELT.data)
  SPELT.data <- remove.young.branches(SPELT.data, age.limit, cut.off)
  SPELT.data <- add.SPELT.contrasts.data(SPELT.data)
  return(SPELT.data)
}

# Fit model of residuals against divergence time for each species pair
fit.lag.model <- function(SPELT.data) {
  lag.model <- lm(SPELT.data$residuals ~ SPELT.data$branch.length)
}

# Overall SPELT function
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

  # Identify cherries (independent pairs of species from one node)
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

# Collate details for summary and plot outputs
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

# Plotting function for SPELT objects
plot.SPELT <- function(SPELT.results) {
  par(bty = "l")
  plot(SPELT.results$data$residuals ~ SPELT.results$data$branch.length, 
       xlab = paste("divergence time (", SPELT.summary.details(SPELT.results)[[2]],")", sep = ""),
       ylab = "residuals", cex.main = 0.8, pch = 16, las = 1,
       main = paste("SPELT results: ", SPELT.summary.details(SPELT.results)[[1]], sep = ""))
  abline(fit.lag.model(SPELT.results$data))
  abline(0,0,lty = 2)
}  

#Summary function for SPELT objects
summary.SPELT <- function(SPELT.results) {
  cat("\nSPELT Details:\n", SPELT.summary.details(SPELT.results)[[1]])
  cat("\n",SPELT.summary.details(SPELT.results)[[2]], "\n")
  print(SPELT.results$summary)
}