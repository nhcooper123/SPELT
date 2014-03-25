# Functions required for SPELT method for detecting evolutionary lag

# Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data) == column.name)
}	
 
# Remove incomplete data
remove.incomplete.data <- function(data, var1.col, var2.col) {
  id <- complete.cases(data[, c(var1.col, var2.col)])
  data <- data[id, ]
  return(data)
}

# ID species in tree that are not in the data
id.missing.tree <- function(phy, data, speciesnames.col) {
  setdiff(phy$tip.label, data[, speciesnames.col])
}

# ID species in data that are not in the tree
id.missing.data <- function(phy, data, speciesnames.col) {
  setdiff(data[, speciesnames.col], phy$tip.label)
}    

# Remove missing species from tree
remove.missing.species.tree <- function(phy, data, speciesnames.col) {
  tree.not.data <- id.missing.tree(phy, data, speciesnames.col)
  if (length(tree.not.data) > 0) {
    phy <- drop.tip(phy, tree.not.data)
  }
  return(phy)
}

# Remove missing species from data
remove.missing.species.data <- function(phy, data, speciesnames.col) {
  data.not.tree <- id.missing.data(phy, data, speciesnames.col)
  if (length(data.not.tree) > 0) {
    matches <- match(data[,speciesnames.col], data.not.tree, nomatch = 0)
    data <- subset(data, matches == 0)
  }
  return(data)
}

# Identify total number of nodes in tree
total.nodes <- function(phy) {
  length(phy$tip.label)
}

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

# Extract data for species
get.data <- function(data, variable.col, species.list) {
  sapply(species.list, function(x) 
         data[which(rownames(data) == x),variable.col])
}

# Extract branch length for contrast between two species
# Need a fix for polytomies!!! == numeric(0)
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
remove.young.branches <- function(SPELT.data, age.limit = NULL) {
  if (!is.null(age.limit)) {
    SPELT.data <- SPELT.data[-(c(which(SPELT.data$branch.length < age.limit))), ]
  }
    if (nrow(SPELT.data) < 3) {
      stop("< 3 branches longer than age limit")
    }
  return(SPELT.data) 
}

# Combining all functions for data collation
get.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, 
                           SPELT.data, age.limit = NULL) {
  SPELT.data <- build.SPELT.data(phy)
  SPELT.data <- add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)
  SPELT.data <- get.raw.contrasts(SPELT.data)
  SPELT.data <- remove.young.branches(SPELT.data, age.limit)
  SPELT.data <- add.SPELT.contrasts.data(SPELT.data)
  return(SPELT.data)
}

# Fit model of residuals against divergence time for each species pair
fit.lag.model <- function(SPELT.data) {
  lag.model <- lm(SPELT.data$residuals ~ SPELT.data$branch.length)
}








# Collate details for summary and plot outputs
SPELT.summary.details <- function(SPELT.results) {
  details <- paste("SPELT results: primary variable = ", SPELT.results$variables$primary.variable,
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
       ylab = "residuals", main = SPELT.summary.details(SPELT.results)[[1]], 
       cex.main = 0.8, pch = 16, las = 1)
  abline(fit.lag.model(SPELT.results$data))
  abline(0,0,lty = 2)
}  

#Summary function for SPELT objects
summary.SPELT <- function(SPELT.results) {
  cat("\nSPELT Details:\n", SPELT.summary.details(SPELT.results)[[1]])
  cat("\n",SPELT.summary.details(SPELT.results)[[2]], "\n")
  print(SPELT.results$summary)
}


