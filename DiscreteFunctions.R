library(ape)
library(geiger) 
library(corHMM)
library(phytools)
library(phangorn)

#You can use code you wrote for the correlation exercise here.

#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?


VisualizeData <- function(phy, data) {  
	print(phy)
	plot(phy)
	print(data)
	names(data)
	str(data)
}


CleanData <- function(phy, data) { treedata(phy,data,sort=TRUE)
}