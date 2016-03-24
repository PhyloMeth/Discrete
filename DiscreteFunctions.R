library(ape)
library(geiger) 
library(corHMM)
library(phytools)
library(phangorn)

#You can use code you wrote for the correlation exercise here.


VisualizeData.discrete <- function(phy, data) {
	#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
	quartz(height=7,width=5)
	layout(matrix(c(1), 1,1)); par(mar=c(0.5,0.5,3,0.5))
	plot.phylo(phy, x.lim=c(0, 2.25), y.lim=c(0,25))
	tiplabels(data $data[,1], frame="n", adj=-15)
	tiplabels(data $data[,2], frame="n", adj=-21)
	text(c(1.7, 2), length(phy$tip.label)+1, c("trait 1", "trait 2"))
}

CleanData <- function(phy, data) {
	treedata(phy, data, sort=TRUE)
}