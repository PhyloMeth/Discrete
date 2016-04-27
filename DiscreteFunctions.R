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

		if (length(unique(c(data[,1],data[,2]))) > 2) {
		quartz()
		plot(data[,1],data[,2])
		traitframe <- as.data.frame(data)
		traitframe[,3] <- rownames(traitframe)
		colnames(traitframe)[3] <- "names"
		quartz()
		color.plot.phylo(phy, traitframe, trait="continuous.data1", taxa.names="names")
		quartz()
		color.plot.phylo(phy, traitframe, trait="continuous.data2", taxa.names="names")
		} else {
		print(table(data[,1], data[,2]))
		quartz()
		barplot(c(table(data[,1]), table(data[,2])))
		quartz()
		plot(phy, label.offset = 3, cex=0.2)
		cols1 <- rep("yellow", length(data[,1]))
		names(cols1) <- row.names(data)
		cols1[data[,1]==1] <- "blue"
		cols2 <- rep("green", length(data[,2]))
		names(cols2) <- row.names(data)
		cols2[data[,2]==1] <- "red"
		tiplabels(pch = 23, col= cols1[phy$tip.label],  bg=cols1[phy$tip.label], adj = 1.4)
		tiplabels(pch = 22, col = cols2[phy$tip.label], bg=cols2[phy$tip.label], adj = 2.5)
		}
}

CleanData <- function(phy, data) {
	#treedata() in Geiger is probably my favorite function in R.
	dataclean <- treedata(phy,data)
	return(dataclean)
}
