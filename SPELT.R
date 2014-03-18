#Function for detecting evolutionary lag
#Using method of Deaner and Nunn 1999

#----FUNCTIONS-----
#Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
}	
 
#Remove incomplete data
remove.incomplete.data <- function(data, var1.col, var2.col) {
  id <- complete.cases(data[,c(var1.col,var2.col)])
  data <- data[id,]
  return(data)
}

#ID species in tree that are not in the data
id.missing.tree <- function(phy, data, speciesnames.col) {
  setdiff(phy$tip.label, data[,speciesnames.col])
}

#ID species in data that are not in the tree
id.missing.data <- function(phy, data, speciesnames.col) {
  setdiff(data[,speciesnames.col], phy$tip.label)
}    

#Remove missing species from tree
remove.missing.species.tree <- function(phy, data, speciesnames.col) {
  tree.not.data <- id.missing.tree(phy, data, speciesnames.col)
  if(length(tree.not.data)>0) {
    phy <- drop.tip(phy, tree.not.data)
  } else {
    phy <- phy
  }
  return(phy)
}

#Remove missing species from data
remove.missing.species.data <- function(phy, data, speciesnames.col) {
  data.not.tree <- id.missing.data(phy, data, speciesnames.col)
  if(length(data.not.tree)>0) {
    matches <- match(data[,speciesnames.col], data.not.tree)
    data <- subset(data, matches!=0)
  } else {
    data <- data
  }
  return(data)
}

#Identify total number of nodes in tree
total.nodes <- function(phy) {
  length(phy$tip.label)
}

#Identify cherries (independent pairs of species from one node)
cherry.nodes <- function(phy) {
  names(which(table(phy$edge[,1][phy$edge[,2]<=total.nodes(phy)])==2))
}

#Exclude branches above an age limit
age.limit <- function(data, branch.col, age.limit)
  data <- data[-(c(which(data[,branch.col]<age.limit))),]
  if(length(data[,branch.col]<3) {
    stop("< 3 branches longer than age.limit")
  } else {
  return(data)
  }
}

#Identify 1st species coming from node
node.species1 <- function(phy, node.list) {
  sapply(node.list, function(x) phy$tip.label[phy$edge[,2][which(phy$edge[,1]==x)][1]])
}

#Identify 2nd species coming from node
node.species2 <- function(phy, node.list) {
  sapply(node.list, function(x) phy$tip.label[phy$edge[,2][which(phy$edge[,1]==x)][2]])
}

#Extract data for species
get.data <- function(data, variable.col, species.list) {
  sapply(species.list, function(x) data[which(rownames(data)==x),variable.col])
}

#Extract branch length for contrast between two species
#Need a fix for polytomies!!! == numeric(0)
branch.length.pair <- function(node.list) {
  sapply(node.list, function(x) phy$edge.length[which(phy$edge[,1]==x)][1])
}

#Build empty dataset for SPELT
build.SPELT.data <- function(phy) {
  SPELT.data <- data.frame(array(dim = c(length(cherry.nodes(phy)),9)))
  names(SPELT.data)<-c("species1", "species2", "species1.primary.var", "species2.primary.var",
        "species1.lag.var","species2.lag.var", "branch.length", "contrast.primary.var", "contrast.lag.var")
}

#Calculate contrasts (primary variable contrast is always positive)

get.raw.contrasts <- function(SPELT.data) { 
  for(i in seq_along(SPELT.data[,1])) {
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

#Fill data set for SPELT
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

#----------------------------
#Actual SPELT function
SPELT<-function(phy, data, primary.variable, lag.variable, speciesnames, age.limit = NA, warn.dropped = TRUE){

  source(SPELT_functions.R)
  require(ape)

------------------------------------------
  #Add checks to code
  #Ensure tree is fully bifurcating
  phy<-multi2di(phy)

  if (!inherits(data, "comparative.data")) 
        stop("data is not a 'comparative' data object.")

#Make rownames into species names

rownames(data)<-data[,colno.speciesnames]

#Warning message showing which species don't match
if(warn.dropped){id.missing.data()

      }
----------------------------------

  #Define variables
  var1.col<-column.ID(data, primary.variable)
  var2.col<-column.ID(data, lag.variable)
  speciesnames.col<-column.ID(data, speciesnames)

  #Tidy up data and tree
  data <- remove.incomplete.data(data, var1.col, var2.col)
  data <- remove.missing.species.data(phy, data, speciesnames.col)
  phy <- remove.missing.species.tree(phy, data, speciesnames.col)

  #Identify cherries (independent pairs of species from one node)
  node.list <- cherry.nodes(phy)
  
  #Build empty dataset for SPELT
  SPELT.data <- build.SPELT.data(phy)
  SPELT.data <- add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)

--------------------
  #Exclude branches above an age limit
  age.limit <- function(data, branch.col, age.limit)
  data <- data[-(c(which(data[,branch.col]<age.limit))),]
  if(length(data[,branch.col]<3) {
    stop("< 3 branches longer than age.limit")
  } else {
  return(data)
  }
}


#-----------------------------------------------------
#Remove branches shorter than user defined age limit
#-----------------------------------------------------

if(!is.na(age.limit)){

		pairs.data<-pairs.data[-(c(which(pairs.data$branch < age.limit))),]
		if(length(pairs.data$branch) < 3){stop("< 3 branches longer than age.limit")}
			}


#--------------------------------------------------

#Fit model of independent contrasts in lag variable
#against independent contrasts in primary variable,
#forced through the origin
#Store residuals in pairs.data

#--------------------------------------------------



model<-lm(contr.lag ~ contr.primary-1, data = pairs.data)#plot contrasts against each other and take residuals

pairs.data$residuals<-unclass(model)$residuals


#--------------------------------------------------

#Fir model of redisuals against divergence time for
#each independent species pair

#--------------------------------------------------


lag.model<-lm(residuals ~ branch, data = pairs.data)



#--------------------------------------------------

#Extract model outputs

#--------------------------------------------------



print(summary(lag.model))



par(bty = "l")

plot(residuals ~ branch, data = pairs.data, xlab = "divergence time", main = "SPELT plot", pch = 16)

abline(lag.model)

abline(0,0,lty = 2)



return(list(data = pairs.data, dropped = list(Tree.not.data = Tree.not.data, Data.not.tree = Data.not.tree), age.limit = age.limit))

}



#--------------------------------------------------

#Example useage
#--------------------------------------------------



library(caper)

data(shorebird)


lag.results<-SPELT(phy = shorebird.tree, data = shorebird.data, primary.variable = F.Mass, lag = Egg.Mass, speciesnames = Species, warn.dropped = TRUE)



str(lag.results)#To look at raw data

#----------------------------------------------------
#With an age limit on the branches of 5 million years
#----------------------------------------------------

lag.results2<-SPELT(phy = shorebird.tree, data = shorebird.data, primary.variable = F.Mass, lag = Egg.Mass, speciesnames = Species, age.limit = 5, warn.dropped = TRUE)



str(lag.results2)