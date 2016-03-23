library(ape)
library(geiger) 
library(corHMM)
library(phytools)
library(phangorn)

#You can use code you wrote for the correlation exercise here.


VisualizeData <- function(data) {
	#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
  pdf("Visualize.pdf")
  par(mfrow=c(1,2))
  plot(data[[1]])
  hist(data[[2]])
  dev.off()
}

CleanData <- function(phy, data) {
	#treedata() in Geiger is probably my favorite function in R.
  treedata(phy,data)
}
