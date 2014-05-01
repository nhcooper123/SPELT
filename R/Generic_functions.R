# Generic functions

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

# Extract data for species
get.data <- function(data, variable.col, speciesnames.col, species.list) {
  sapply(species.list, function(x) 
         data[which(data[, speciesnames.col] == x),variable.col])
}