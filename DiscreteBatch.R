#You can use code you wrote for the correlation exercise here.
setwd("~/Desktop/PhyloMeth/DiscreteTraits")

source("DiscreteFunctions.R")

tree <- read.tree("ZanneEtAlChronogram1.tre")

dat<-c("Daucus_carota" , "Conyza_canadensis" , "Setaria_viridis" , "Andropogon_gerardii" , "Pueraria_montana" , "Ambrosia_trifida" , "Datura_stramonium" , "Smallanthus_uvedalia" , "Perilla_frutescens" , "Verbesina_jacksonii" , "Vernonia_gigantea" , "Lespedeza_cuneata" , "Ambrosia_artemisiifolia" , "Cichorium_intybus" , "Albizia_julibrissin" , "Cercis_canadensis" , "Sorghum_halepense" , "Sorghastrum_nutans" , "Ligustrum_sinense" , "Gleditsia_triacanthos" , "Microstegium_vimineum" , "Cyperus_strigosus" , "Mosla_chinensis" , "Symphyotrichum_pilosum" , "Erechtites_hieraciifolius" , "Lespedeza_bicolor" , "Helenium_autumnale" , "Trifolium_pratense" , "Acer_saccharum" , "Lonicera_maackii" , "Liriodendron_tulipifera" , "Ailanthus_altissima")

#discrete.data <- read.csv(file="____PATH_TO_DATA_OR_SOME_OTHER_WAY_OF_GETTING_TRAITS____", stringsAsFactors=FALSE) #death to factors.

invasive<-c(1,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1)
woody<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,1,1,1,1)

#discrete.data<-cbind.data.frame(invasive,woody)
#row.names(discrete.data)<-dat

names(invasive)<-dat
names(woody)<-dat

cleaned.discrete.I <- CleanData(tree, invasive)
cleaned.discrete.W <- CleanData(tree, woody)

VisualizeData(cleaned.discrete.I$phy, invasive)
VisualizeData(cleaned.discrete.W$phy, woody)

#First, let's use parsimony to look at ancestral states
invasive.phyDat <- phyDat(invasive, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn
anc.pI <- ancestral.pars(cleaned.discrete.I$phy, invasive.phyDat)
plotAnc(cleaned.discrete.I$phy, anc.pI, 1)

woody.phyDat <- phyDat(woody, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn
anc.pW <- ancestral.pars(cleaned.discrete.W$phy, woody.phyDat)
plotAnc(cleaned.discrete.W$phy, anc.pW, 1)

#Do you see any uncertainty? What does that meean for parsimony?

#now plot the likelihood reconstruction
anc.ml.I <- ancestral.pml(pml(cleaned.discrete.I$phy, invasive.phyDat), type="ml")
plotAnc(cleaned.discrete.I$phy, anc.ml.I, 1)

anc.ml.W <- ancestral.pml(pml(cleaned.discrete.W$phy, woody.phyDat), type="ml")
plotAnc(cleaned.discrete.W$phy, anc.ml.W, 1)
#How does this differ from parsimony? 
#Why does it differ from parsimony?
#What does uncertainty mean?

#How many changes are there in your trait under parsimony? 
parsimony.score.I <- parsimony(cleaned.discrete.I$phy, invasive.phyDat)
print(parsimony.score.I)

parsimony.score.W <- parsimony(cleaned.discrete.W$phy, woody.phyDat)
print(parsimony.score.W)
#Can you estimate the number of changes under a likelihood-based model? 

#Well, we could look at branches where the reconstructed state changed from one end to the other. But that's not really a great approach: at best, it will underestimate the number of changes (we could have a change on a branch, then a change back, for example). A better approach is to use stochastic character mapping.

estimated.histories.I <- make.simmap(cleaned.discrete.I$phy, invasive, model="ARD", nsim=5)
estimated.histories.W <- make.simmap(cleaned.discrete.W$phy, woody, model="ARD", nsim=5)

#always look to see if it seems reasonable
plotSimmap(estimated.histories.I)
plotSimmap(estimated.histories.W)

counts.I <- countSimmap(estimated.histories.I)
print(counts.I)

counts.W <- countSimmap(estimated.histories.W)
print(counts.W)

#Depending on your biological question, investigate additional approaches:
#  As in the correlation week, where hypotheses were examined by constraining rate matrices, one can constrain rates to examine hypotheses. corHMM, ape, and other packages have ways to address this.
#  Rates change over time, and this could be relevant to a biological question: have rates sped up post KT, for example. Look at the models in geiger for ways to do this.
#  You might observe rates for one trait but it could be affected by some other trait: you only evolve wings once on land, for example. corHMM can help investigate this.