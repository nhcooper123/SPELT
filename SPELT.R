#SPELT method for detecting evolutionary lag
#Based on Deaner and Nunn 1999

SPELT <- function(phy, data, primary.variable, lag.variable, 
                  speciesnames, age.limit = NA, warn.dropped = TRUE) {

  source(SPELT_functions.R)
  require(ape)

#------------------------------------------
  #Add checks to code. Check phy = phylo object etc.
# if (!inherits(data, "comparative.data")) 
#        stop("data is not a 'comparative' data object.")
#------------------------------------------------------
  #Define variables
  var1.col<-column.ID(data, primary.variable)
  var2.col<-column.ID(data, lag.variable)
  speciesnames.col<-column.ID(data, speciesnames)

  #Tidy up data and tree
  data <- remove.incomplete.data(data, var1.col, var2.col)
  phy <- multi2di(phy)
  data <- remove.missing.species.data(phy, data, speciesnames.col)
  phy <- remove.missing.species.tree(phy, data, speciesnames.col)

  if(warn.dropped){
    tree.not.data <- id.missing.tree(phy, data, speciesnames.col)
    data.not.tree <- id.missing.data(phy, data, speciesnames.col)
      }

  #Identify cherries (independent pairs of species from one node)
  node.list <- cherry.nodes(phy)
  
  #Build empty dataset for SPELT
  SPELT.data <- build.SPELT.data(phy)
  SPELT.data <- add.SPELT.data(phy, data, node.list, var1.col, var2.col, SPELT.data)

  #Remove branches shorter than user defined age limit
  if(!is.na(age.limit)) {
    SPELT.data<-SPELT.data[-(c(which(SPELT.data[,7]<age.limit))),]
		if(length(pairs.data$branch) < 3) {
      stop("<3 branches longer than age limit")
    }		
  }

  #Fit models and extract residuals
  SPELT.data[,10] <- contrasts.model.residuals(SPELT.data)

  #Fit SPELT model 
  SPELT.model <- lag.model(SPELT.data)
  
  #Outputs
  SPELT.results <- list(call = SPELT.model$call, variables = list(primary.variable = primary.variable, 
                        lag.variable = lag.variable), data = SPELT.data, age.limit = age.limit,
                        dropped = list(tree.not.data = tree.not.data, data.not.tree = data.not.tree), 
                        Nnodes = length(SPELT.data[,1])) 

  class(SPELT.results) <- "SPELT"

  return(SPELT.results)
}