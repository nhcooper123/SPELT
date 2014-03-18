#Examples for using SPELT with shorebird data from Thomas et al. 2006
#Does egg mass evolution lag behind female mass?

library(caper)
data(shorebird)

SPELT.result <- SPELT(phy = shorebird.tree, data = shorebird.data, 
	                 primary.variable = F.Mass, lag = Egg.Mass, 
	                 speciesnames = Species, warn.dropped = TRUE)

summary(SPELT.result)
print(SPELT.result)

#With an age limit on the branches of 5 million years

SPELT.result.5MY <- SPELT(phy = shorebird.tree, data = shorebird.data, 
	                 primary.variable = F.Mass, lag = Egg.Mass, 
	                 speciesnames = Species, age.limit = 5, warn.dropped = TRUE)

summary(SPELT.result.5MY)
print(SPELT.result.MY)