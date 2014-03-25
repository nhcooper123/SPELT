#--------------------------------------------------
#SPELT method
#--------------------------------------------------

#Requires ape


#--------------------------------------------------

#Input
#--------------------------------------------------


#primary.variable = primary variable which you think is driving evolution of the lagged trait

#lag.variable = variable you think is lagging behind the primary variable

#phy = phylogeny

#data = dataset as a data.frame

#speciesnames = name of the column in your dataset with the species names

#age.limit = numeric variable allowing the user to set a minimum age for branches to use
#This will result in an error with < 3 branches of sufficient length because it is not 
	#possible to fit the model with less than 3 species pairs
#warn.dropped = if TRUE gives a list of species dropped because they are not in the tree and data

#--------------------------------------------------

#Output
#--------------------------------------------------

#1 - plot of residuals versus divergence time

#2 - model output of residuals versus divergence time

#3 - dataset with summary variables

	#i.e. pair number, the species in each pair, values of their response and lag variables, 

	#branch length, contrasts and residuals from model
	#Also $dropped showing the species in the tree and not the data
	#or in the data and not the tree.


#--------------------------------------------------

#SPELT function
#--------------------------------------------------

SPELT<-function(phy, data, primary.variable, lag.variable, speciesnames, age.limit = NA, warn.dropped = TRUE){



require(ape)


#--------------------------------------------------

#Modify variables

#--------------------------------------------------


primary.variable<-deparse(substitute(primary.variable))

lag.variable<-deparse(substitute(lag.variable))

speciesnames<-deparse(substitute(speciesnames))



#--------------------------------------------------

#Identify columns for variables

#--------------------------------------------------



colno.primary.variable<-which(names(data)==primary.variable)

colno.lag<-which(names(data)==lag.variable)

colno.speciesnames<-which(names(data)==speciesnames)


#-------------------------------------------------
#Warning message showing which species don't match
#-------------------------------------------------

if(warn.dropped){

   Tree.not.data <- setdiff(phy$tip.label, data[,colno.speciesnames])
   Data.not.tree <- setdiff(data[,colno.speciesnames], phy$tip.label)
		
			}


#--------------------------------------------------

#Strip data and tree to remove missing values

#--------------------------------------------------



id<-complete.cases(data[, c(colno.primary.variable,colno.lag)])

data<-data[id,]



matches<-match(phy$tip.label, data[,colno.speciesnames], nomatch = 0)  

not<-subset(phy$tip.label, matches == 0)



if(length(not)>0){

phy<-drop.tip(phy, not)

}else{phy<-phy}



matches2<-match(data[,colno.speciesnames], phy$tip.label, nomatch = 0)  

data<-subset(data, matches2 !=0)

#--------------------------------------------------

#Ensure tree is fully bifurcating

#--------------------------------------------------

phy<-multi2di(phy)


#--------------------------------------------------

#Make rownames into species names

#--------------------------------------------------



rownames(data)<-data[,colno.speciesnames]



#--------------------------------------------------

#Identify independent species pairs

#--------------------------------------------------



total.nodes <- length(phy$tip.label)

cherry.nodes<-names(which(table(phy$edge[, 1][phy$edge[, 2] <= total.nodes])==2)) #detects independent pairs (cherries)



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