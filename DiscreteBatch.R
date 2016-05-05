#You can use code you wrote for the correlation exercise here.
setwd("~/GitHub/DiscreteTraits")
source("DiscreteFunctions.R")


tree<-pbtree(n=20,scale=1,tip.label=(c("T.ochroleucum","T.palmeri","T.calocephalum","T.nanum","T.pallescens","T.repens","T.amabile","T.decorum","T.wigginsii","T.virginicum","T.vesiculosum","T.vernum","T.uniflorum","T.tomentosa","T.thalii","T.sylvatica","T.alpestre","T.argentinense","T.pratense","T.hybridum")))

dis.1<-c(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1)
dis.2<-c(0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1)
discrete.data<-cbind.data.frame(row.names=tree$tip.label,dis.1,dis.2)

# simulating data instead
#tree <- read.tree("____PATH_TO_TREE_OR_SOME_OTHER_WAY_OF_GETTING_A_TREE____")
#discrete.data <- read.csv(file="____PATH_TO_DATA_OR_SOME_OTHER_WAY_OF_GETTING_TRAITS____", stringsAsFactors=FALSE) #death to factors.

cleaned.discrete <- CleanData(tree, discrete.data)
VisualizeData(tree, cleaned.discrete)


#First, let's use parsimony to look at ancestral states
cleaned.discrete.phyDat <- phyDat(cleaned.discrete$data, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn
anc.p <- ancestral.pars(cleaned.discrete$phy, cleaned.discrete.phyDat)
plotAnc(cleaned.discrete$phy, anc.p, 1)

#Do you see any uncertainty? What does that mean for parsimony?
#yes some states are uncertain.  This means that either state is equally parsimonious.


#now plot the likelihood reconstruction
anc.ml <- ancestral.pml(pml(cleaned.discrete$phy, cleaned.discrete.phyDat), type="ml")
plotAnc(tree, anc.ml, 1)

#How does this differ from parsimony? 
# The node states are not 50/50 if uncertain
# they give a percentage chance for each state
#Why does it differ from parsimony?
# it plots the most likely state
#What does uncertainty mean?
# that either state has a given percent chance of being correct at a given node

#How many changes are there in your trait under parsimony? 
parsimony.score <- parsimony(cleaned.discrete$phy, cleaned.discrete.phyDat)
print(parsimony.score)

#Can you estimate the number of changes under a likelihood-based model? 
# not really
#Well, we could look at branches where the reconstructed state changed from one end to the other. But that's not really a great approach: at best, it will underestimate the number of changes (we could have a change on a branch, then a change back, for example). A better approach is to use stochastic character mapping.

estimated.histories <- make.simmap(cleaned.discrete$phy, cleaned.discrete$data[,1], model="ARD", nsim=5)

#always look to see if it seems reasonable
plotSimmap(estimated.histories)

counts <- countSimmap(estimated.histories)
print(counts)

#Depending on your biological question, investigate additional approaches:
#  As in the correlation week, where hypotheses were examined by constraining rate matrices, one can constrain rates to examine hypotheses. corHMM, ape, and other packages have ways to address this.
#  Rates change over time, and this could be relevant to a biological question: have rates sped up post KT, for example. Look at the models in geiger for ways to do this.
#  You might observe rates for one trait but it could be affected by some other trait: you only evolve wings once on land, for example. corHMM can help investigate this.