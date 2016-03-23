#You can use code you wrote for the correlation exercise here.
source("DiscreteFunctions.R")
tree <- read.tree("Eurycea_Tree")

discrete.1<-round(runif(length(tree$tip.label)))
discrete.2<-round(runif(length(tree$tip.label)))
names(discrete.1) <- tree$tip.label
names(discrete.2) <- tree$tip.label
# discrete.data <- data.frame(discrete.1,discrete.2) # In case I later want to put these data in a dataframe together.
# row.names(discrete.data) <- tree$tip.label # In case I later want to put these data in a dataframe together.

cleaned.discrete <- CleanData(tree, discrete.1)

VisualizeData(cleaned.discrete)

#First, let's use parsimony to look at ancestral states
cleaned.discrete.phyDat <- phyDat(cleaned.discrete$data, type="USER",levels=c(0,1)) #phyDat is a data format used by phangorn
anc.p <- ancestral.pars(tree, cleaned.discrete.phyDat)
plotAnc(tree, anc.p, 1)

#Do you see any uncertainty? What does that meean for parsimony?
# There is uncertainty at some of the nodes. It means that there are multiple, equally 'likely'
# (i.e., smallest number of changes) ways to reconstruct these states.

#now plot the likelihood reconstruction
anc.ml <- ancestral.pml(pml(tree, cleaned.discrete.phyDat), type="ml")
plotAnc(tree, anc.ml, 1)

#How does this differ from parsimony? 
# The nodal estimates changed, and they're no longer 50/50 when uncertainty exists.
#Why does it differ from parsimony?
# The model is different; it incorporates rates of change along branchs of varying lengths.
#What does uncertainty mean?
# It represents the probability of those states.

#How many changes are there in your trait under parsimony? 
parsimony.score <- parsimony(tree, cleaned.discrete.phyDat)
print(parsimony.score)

#Can you estimate the number of changes under a likelihood-based model? 

#Well, we could look at branches where the reconstructed state changed from one end to the other. But that's not really a great approach: at best, it will underestimate the number of changes (we could have a change on a branch, then a change back, for example). A better approach is to use stochastic character mapping.

estimated.histories <- make.simmap(tree, discrete.1, model="ARD", nsim=5)

#always look to see if it seems reasonable
plotSimmap(estimated.histories)

counts <- countSimmap(estimated.histories)
print(counts)

#Depending on your biological question, investigate additional approaches:
#  As in the correlation week, where hypotheses were examined by constraining rate matrices, one can constrain rates to examine hypotheses. corHMM, ape, and other packages have ways to address this.
#  Rates change over time, and this could be relevant to a biological question: have rates sped up post KT, for example. Look at the models in geiger for ways to do this.
#  You might observe rates for one trait but it could be affected by some other trait: you only evolve wings once on land, for example. corHMM can help investigate this.