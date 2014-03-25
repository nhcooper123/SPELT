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
  names(which(table(phy$edge[, 1][phy$edge[, 2]<=total.nodes(phy)]) == 2))
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
  sapply(node.list, function(x) phy$edge.length[which(phy$edge[,1]==x)][1])
}

# Build empty dataset for SPELT
build.SPELT.data <- function(phy) {
  SPELT.data <- data.frame(array(dim = c(length(cherry.nodes(phy)),10)))
  colnames(SPELT.data)<-c("species1", "species2", "species1.var1", 
                      "species2.var1", "species1.var2",
                       "species2.var2", "branch.length", 
                      "contrast.var1", "contrast.var2", 
                       "residuals")
}

# Calculate contrasts (primary variable contrast is always positive)
get.raw.contrasts <- function(SPELT.data) { 
  for(i in seq_along(SPELT.data$)) {
    if(SPELT.data[i,3] > SPELT.data[i,4]) {
      SPELT.data[i,8] <- SPELT.data[i,3] - SPELT.data[i,4]
      SPELT.data[i,9] <- SPELT.data[i,5] - SPELT.data[i,6]
    } else {
      SPELT.data[i,8] <- SPELT.data[i,4] - SPELT.data[i,3]
      SPELT.data[i,9] <- SPELT.data[i,6] - SPELT.data[i,5]
    }
  }  
  return(SPELT.data)
}

# Fill data set for SPELT
add.SPELT.data <- function(phy, data, node.list, var1.col, var2.col, SPELT.data) {
  SPELT.data[,1] <- node.species1(phy,node.list)
  SPELT.data[,2] <- node.species2(phy,node.list)
  SPELT.data[,3] <- get.data(data, var1.col, SPELT.data[,1])
  SPELT.data[,4] <- get.data(data, var1.col, SPELT.data[,2])
  SPELT.data[,5] <- get.data(data, var2.col, SPELT.data[,1])
  SPELT.data[,6] <- get.data(data, var2.col, SPELT.data[,2])
  SPELT.data[,7] <- branch.length.pair(node.list)
  get.raw.contrasts(SPELT.data)
}

# Fit model of lag variable contrasts ~ primary variable contrasts (forced through the origin)
contrasts.model.residuals <- function(SPELT.data) {
  model<-lm(SPELT.data[,9] ~ SPELT.data[,8] - 1)
  residuals<-unclass(model)$residuals
  return(residuals)
}

#Fit model of redisuals against divergence time for each species pair
lag.model <- function(SPELT.data) {
  lag.model <- lm(SPELT.data[,10] ~ SPELT.data[,7])
}

#Plotting function for SPELT objects
plot.SPELT <- function(SPELT.results) {
  plot(SPELT.data[,10] ~ SPELT.data[,7], xlab = "divergence time", 
       ylab = "residuals", main = "SPELT plot", pch = 16)
  abline(lag.model(SPELT.data))
  abline(0,0,lty = 2)
}  

#Summary function for SPELT objects
summary.SPELT <- function(object, ...) {
  ans <- list(call = object$call)
  class(ans) <- "summary.SPELT"
  estimate <- unclass(summary(object)$coefficients)[1:2]
  sterr <- unclass(summary(object)$coefficients)[3:4]
  t <- unclass(summary(object)$coefficients)[5:6]
  p <- unclass(summary(object)$coefficients)[7:8]
  coef <- cbind(estimate, sterr, t, p)
  colnames(coef) <- c("Estimate", "Std. Error", "t value", 
                      "Pr(>|t|)")
  ans$coefficients <- coef
  ans$df <- unclass(summary(object)$df)[2]
  ans$AIC <- AIC(object)
  ans$r.squared <- unclass(summary(object)$r.squared)
  return(ans)
  ans$fitted <- fitted(object)
  ans$residuals <- residuals(object)
  
}