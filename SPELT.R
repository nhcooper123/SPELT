#Function for detecting evolutionary lag
#Using method of Deaner and Nunn 1999

#----FUNCTIONS-----
#Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
}	

#Check species match
id.missing.data <- function(phy, data, speciesnames) {
  speciesnames.col <- column.ID(data, speciesnames)
  tree.not.data <- setdiff(phy$tip.label, data[,speciesnames.col])
  data.not.tree <- setdiff(data[,speciesnames.col], phy$tip.label)
  return(tree.not.data)
  return(data.not.tree)
}  		

#Remove incomplete data
remove.incomplete.data <- function(data, var1, var2) {
  id<-complete.cases(data[,c(var1,var2)])
  data<-data[id,]
  return(data)
}

#Remove missing species
remove.missing.data <- function(phy, data, speciesnames) {
  id.missing.data(phy, data, speciesnames)
  if(length(tree.not.data)>0) {
    phy <- drop.tip(phy, tree.not.data)
  } 
  else {
    phy <- phy
  }
  if(length(data.not.tree)>0) {
    speciesnames.col <- column.ID(data, speciesnames)
    matches<-match(data[,speciesnames.col], data.not.tree)
    data <- subset(data, matches!=0)
  } 
  else {
    data <- data
  }
  return(data)
  return(phy)
}

#Identify total number of nodes in tree
total.nodes <- function(phy) {
  length(phy$tip.label)
}

#Exclude branches above an age limit
age.limit <- function(data, branch.col, age.limit)
  branch <- column.ID(data, branch.col)
  }
  data <- data[-(c(which(data[,branch] < age.limit))),]
  if(length(data[,branch]<3) {
    stop("< 3 branches longer than age.limit")
}

#Identify cherries (independent pairs of nodes)
cherry.nodes <- function(phy) {
  names(which(table(phy$edge[,1][phy$edge[,2] <= total.nodes])==2))
}

#Bag
for (i in 1:length(cherry.nodes)){

	pairs.data$pair[i]<-i	

	pairs.data$species1[i]<-phy$tip.label[phy$edge[,2][which(phy$edge[,1]==cherry.nodes[i])][1]]#species 1 of the pair

	pairs.data$species2[i]<-phy$tip.label[phy$edge[,2][which(phy$edge[,1]==cherry.nodes[i])][2]]#species 2 of the pair

	pairs.data$species1_primary[i]<-data[which(rownames(data)==pairs.data$species1[i]),colno.primary.variable]#primary variable for species 1

	pairs.data$species2_primary[i]<-data[which(rownames(data)==pairs.data$species2[i]),colno.primary.variable]#primary variable for species 2

	pairs.data$species1_lag[i]<-data[which(rownames(data)==pairs.data$species1[i]),colno.lag]#lag variable for species 1

	pairs.data$species2_lag[i]<-data[which(rownames(data)==pairs.data$species2[i]),colno.lag]#lag variable for species 2

	pairs.data$branch[i]<-phy$edge.length[which(phy$edge[,1]==cherry.nodes[i])][1]#branch length for contrast




#----------------------------
#Actual SPELT function
SPELT<-function(phy, data, primary.variable, lag.variable, speciesnames, age.limit = NA, warn.dropped = TRUE){

source(SPELT_functions.R)
require(ape)

#Add checks to code
#Ensure tree is fully bifurcating
phy<-multi2di(phy)

if (!inherits(data, "comparative.data")) 
        stop("data is not a 'comparative' data object.")

#Define variables
primary.var.col<-column.ID(data, primary.variable)
lag.var.col<-column.ID(data, lag.variable)
speciesnames.col<-column.ID(data, speciesnames)

#Warning message showing which species don't match
if(warn.dropped){id.missing.data()

			}

#Strip data and tree to remove missing values


#Make rownames into species names

rownames(data)<-data[,colno.speciesnames]


#Identify independent species pairs

#--------------------------------------------------




#--------------------------------------------------

#Setup dataframe for results

#--------------------------------------------------



pairs.data<-data.frame(array(dim = c(length(cherry.nodes),10)))#makes new dataframe to put summary data into

names(pairs.data)<-c("pair", "species1", "species2", "species1_primary", "species2_primary",
				"species1_lag","species2_lag", "branch", "contr.primary", "contr.lag")



#---------------------------------------------------------------------------------

#Identify species, variables and branch lengths for each independent species pair

#Input into pairs.data dataframe

#---------------------------------------------------------------------------------



for (i in 1:length(cherry.nodes)){



	pairs.data$pair[i]<-i	

	pairs.data$species1[i]<-phy$tip.label[phy$edge[,2][which(phy$edge[,1]==cherry.nodes[i])][1]]#species 1 of the pair

	pairs.data$species2[i]<-phy$tip.label[phy$edge[,2][which(phy$edge[,1]==cherry.nodes[i])][2]]#species 2 of the pair

	pairs.data$species1_primary[i]<-data[which(rownames(data)==pairs.data$species1[i]),colno.primary.variable]#primary variable for species 1

	pairs.data$species2_primary[i]<-data[which(rownames(data)==pairs.data$species2[i]),colno.primary.variable]#primary variable for species 2

	pairs.data$species1_lag[i]<-data[which(rownames(data)==pairs.data$species1[i]),colno.lag]#lag variable for species 1

	pairs.data$species2_lag[i]<-data[which(rownames(data)==pairs.data$species2[i]),colno.lag]#lag variable for species 2

	pairs.data$branch[i]<-phy$edge.length[which(phy$edge[,1]==cherry.nodes[i])][1]#branch length for contrast

				

#--------------------------------------------------							

#Calculate independent contrasts, maintaining sign

#Primary variable contrast is always positive

#Lag variable contrast can be positive or negative

#--------------------------------------------------


if(pairs.data$species1_primary[i] > pairs.data$species2_primary[i]){

		pairs.data$contr.primary[i]<-	pairs.data$species1_primary[i] - pairs.data$species2_primary[i]						

		pairs.data$contr.lag[i]<-pairs.data$species1_lag[i] - pairs.data$species2_lag[i]		

		}else{	

		pairs.data$contr.primary[i]<-	pairs.data$species2_primary[i] - pairs.data$species1_primary[i]

		pairs.data$contr.lag[i]<-pairs.data$species2_lag[i] - pairs.data$species1_lag[i]

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