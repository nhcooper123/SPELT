library(devtools)
install_github("SPELT", username = "nhcooper123")
library(SPELT)

#To look at the helpfiles
help(package = "SPELT") #This has a whole bunch of functions that are internal to SPELT
?SPELT #This will show you the help for the main function

datamain <- read.csv("/Users/charlienunn/Documents/Research/Current Projects/Evolutionary Lag COOPER/Empirical - Number of males/Lindenfors_et_al_2004_Data_FINAL.csv")
head(datamain)

primateTreeBlock <- read.nexus("~/Documents/Research/Current Projects/Evolutionary Lag COOPER/Empirical - Number of males/consensusTree_10kTrees_Primates_Version3.ADD3TAXA.nex")
primateTree <- primateTreeBlock[[1]]
plot(primateTree)

#Running SPELT
lag.results <- SPELT(primateTree, datamain, "LogFemales", "LogMales", "Species")

#Looking at the results:

#This gives you the model summary you need
summary(lag.results)
#This plots the results
plot(lag.results)
#You can look at the raw contrasts data here
str(lag.results$data)
#The full output can be seen using this
lag.results

#If you want mins and maxes and means:
min(lag.results$data$branch.length)
max(lag.results$data$branch.length)
mean(lag.results$data$branch.length)

#There's no need to do the model manually. SPELT does it for you. And also provides the summary
#but if you wanted to do this...
summary(glm(lag.results$data$contrast.var2 ~ lag.results$data$contrast.var1 -1))

#And with branch length limits:
lag.results.10MY <- SPELT(primateTree, datamain, "LogFemales", "LogMales", "Species", age.limit = 10)
summary(lag.results.10MY)

lag.results.5MY <- SPELT(primateTree, datamain, "LogFemales", "LogMales", "Species", age.limit = 5)
summary(lag.results.5MY)


