#You can use code you wrote for the correlation exercise here.
source("DiscreteFunctions.R")
tree <- read.tree("eucTree.tre")
tree$tip.label <- gsub("_"," ", tree$tip.label); tree$tip.label[8] <- "E. tenuiramis"
treeTransform <- tree; treeTransform$edge.length <- (tree$edge.length)^(1/4) # transform, just to be able to see speciation events


discrete.data <- read.csv(file="eucDataDiscretized.csv", stringsAsFactors=FALSE, row.names=1) #death to factors.

cleaned.discrete <- CleanData(tree, discrete.data)

VisualizeData.discrete(tree, cleaned.discrete)

#First, let's use parsimony to look at ancestral states
cleaned.discrete.phyDat <- phyDat(cleaned.discrete$data, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn
anc.p <- ancestral.pars(tree, cleaned.discrete.phyDat)
plotAnc(treeTransform, anc.p, 1)

#Do you see any uncertainty? What does that meean for parsimony?
	#no uncertainty because this is the tree with the least number of character state transitions

#now plot the likelihood reconstruction
anc.ml <- ancestral.pml(pml(tree, cleaned.discrete.phyDat), type="ml")
plotAnc(treeTransform, anc.ml,1)

#How does this differ from parsimony? 
#Why does it differ from parsimony?
#What does uncertainty mean?

# This is different from parsimony because it doesn't pick a single set of node states, it uses likelihood to determine a set of probabilities for states at each node. Likelihood is the probability of the data given the reconstruction; multiple reconstructions can have some likelihood meaning that there can be uncertainty

#How many changes are there in your trait under parsimony? 
parsimony.score <- parsimony(tree, cleaned.discrete.phyDat)
print(parsimony.score)

#Can you estimate the number of changes under a likelihood-based model? 


#Well, we could look at branches where the reconstructed state changed from one end to the other. But that's not really a great approach: at best, it will underestimate the number of changes (we could have a change on a branch, then a change back, for example). A better approach is to use stochastic character mapping.

estimated.histories <- make.simmap(tree, cleaned.discrete$data[,1], model="ARD", nsim=5)

#always look to see if it seems reasonable
plotSimmap(estimated.histories)

counts <- countSimmap(estimated.histories)
print(counts)

#Depending on your biological question, investigate additional approaches:
#  As in the correlation week, where hypotheses were examined by constraining rate matrices, one can constrain rates to examine hypotheses. corHMM, ape, and other packages have ways to address this.
#  Rates change over time, and this could be relevant to a biological question: have rates sped up post KT, for example. Look at the models in geiger for ways to do this.
#  You might observe rates for one trait but it could be affected by some other trait: you only evolve wings once on land, for example. corHMM can help investigate this.





# Are traits correlated?

cleaned.discrete$dataCor <- data.frame (species=cleaned.discrete$phy$tip.label, cleaned.discrete$data)

redER <- corDISC(cleaned.discrete$phy, cleaned.discrete$dataCor, ntraits=2, rate.mat=nocor, model="ER")
redARD <- corDISC(cleaned.discrete$phy, cleaned.discrete$dataCor, ntraits=2, rate.mat=cor, model="ER")

redER $loglik; redARD $loglik # models assuming equal rates (traits are independent of one another) and different rates (traits are correlated) are no different
