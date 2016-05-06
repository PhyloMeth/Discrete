library(ape)
library(geiger)
library(corHMM)
library(phytools)
library(phangorn)

#You can use code you wrote for the correlation exercise here.


VisualizeData <- function(phy, data) {
	#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
is.rooted(phy)
is.binary.tree(phy)
is.ultrametric(phy)

print(data)

	if (length(unique(c(data[,1]))) > 2) {
	quartz()
	plot(data[,1])
	traitframe <- as.data.frame(data)
	traitframe[,2] <- rownames(traitframe)
	colnames(traitframe)[1] <- "continuous.data"
	colnames(traitframe)[2] <- "names"
	quartz()
	color.plot.phylo(phy, traitframe, trait="continuous.data", taxa.names="names")
	} else {
	print(table(data[,1]))
	quartz()
	barplot(c(table(data[,1])))
	quartz()
	plot(phy, label.offset = 3, cex=0.2)
	cols1 <- rep("yellow", length(data[,1]))
	names(cols1) <- row.names(data)
	cols1[data[,1]==1] <- "blue"
	tiplabels(pch = 23, col= cols1[phy$tip.label],  bg=cols1[phy$tip.label], adj = 1.4)
	}
}

CleanData <- function(phy, data) {
	#treedata() in Geiger is probably my favorite function in R.
	dataclean <- treedata(phy,data)
	return(dataclean)
}

RunDeltaModel<-function(timetransform, phy, data) {
	print(paste("Now starting model",timetransform))
	return(fitDiscrete(phy, data, model="ARD", transform=timetransform))
}

RunCorrModel<-function(correlated, phy, data) {
	print(paste("Now starting model",correlated))
	return(corDISC(phy, data, ntraits=2, model=correlated, node.states="joint"))
}
