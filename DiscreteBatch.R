#You can use code you wrote for the correlation exercise here.
source("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/DiscreteTraits/DiscreteFunctions.R")
tree <- read.tree("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/Correlation/dateddung2.tre")
#discrete.data <- read.csv(file="____PATH_TO_DATA_OR_SOME_OTHER_WAY_OF_GETTING_TRAITS____", stringsAsFactors=FALSE) #death to factors.


# discrete.data1 <- rTraitDisc(tree, model = "ER", k=2, states=c(0,1))
# discrete.data2 <- rTraitDisc(tree, model = "ARD", k=2, rate=c(0.3, 0.5), states=c(0,1))
# discrete.data <- cbind(discrete.data1, discrete.data2)
# discrete.data[discrete.data==1] <- 0
# discrete.data[discrete.data==2] <- 1
# discrete.data1[discrete.data1==1] <- 0
# discrete.data1[discrete.data1==2] <- 1

discrete.data <- rTraitDisc(tree, model = "ER", k=2, states=c(0,1))


cleaned.discrete <- CleanData(tree, discrete.data)

VisualizeData(tree, cleaned.discrete$data)

#First, let's use parsimony to look at ancestral states
cleaned.discrete.phyDat <- phyDat(cleaned.discrete$data, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn
anc.p <- ancestral.pars(tree, cleaned.discrete.phyDat)
quartz()
plotAnc(tree, anc.p, 1, cex.pie=0.3)

#Do you see any uncertainty? What does that meean for parsimony?

#Yes, but only few (2 nodes) have half-half filled pie charts, meaning they could be either state. Parsimony goes for minimal number of steps, so it will decide on unambiguous states wherever it can.

#now plot the likelihood reconstruction
anc.ml <- ancestral.pml(pml(tree, cleaned.discrete.phyDat), type="ml")
quartz()
plotAnc(tree, anc.ml, 1, cex.pie=0.3)

#How does this differ from parsimony?
# There is a lot more uncertainty, basically all the deeper nodes are 50:50 ambiguous and a trait only gets higher once we move closer to clades with homogenious tip traits.
#Why does it differ from parsimony?
# It estimates the likelihood of a trait state based on a model of trait evolution, making room for different possible results which may have a different likelihood.
#What does uncertainty mean?
# That we cannot be entirely sure of the result we got, but that some results are more likely than others

#How many changes are there in your trait under parsimony?
parsimony.score <- parsimony(tree, cleaned.discrete.phyDat)
print(parsimony.score)

#Can you estimate the number of changes under a likelihood-based model?
#(are you answering the question yourself below?) You could define a cutoff value (e.g. 0.5) and see how often it crosses that, but that wouldn't necessarily be realistic...

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
