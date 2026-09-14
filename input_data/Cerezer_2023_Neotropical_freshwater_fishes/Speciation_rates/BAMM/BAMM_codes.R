#### :::: Codes used to calculate and extract tip speciation rates from BAMM :::: ####


# Required packages --------------------------------------------------------
library(ape)
library(BAMMtools)
library(ade4)
library(magrittr)
library(coda)
library(pander)
library(utils)
library(phylobase)
library(ggplot2)
library(phytools)


# Importing phylogeny -----------------------------------------------------
tree <- read.tree("../../Datasets/Phylogeny/NFF_subset_2638spp.tre")

plot(tree,cex=0.05)
axisPhylo()


# Set priors --------------------------------------------------------------
PreviousPriors <- setBAMMpriors(phy = tree)


# Assessing MCMC convergence ----------------------------------------------
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.20 * nrow(mcmcout)) #discard the first 20% of samples as burnin
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

plot(postburn$logLi~postburn$generation,col="red")
plot(postburn$N_shifts~postburn$generation,col="blue")

effectiveSize(postburn$N_shifts) # Effective sample size should ideally be higher than 200
effectiveSize(postburn$logLik) # Effective sample size should ideally be higher than 200


# Importing BAMM outputs --------------------------------------------------
ed<- getEventData(tree, "event_data.txt",burnin = 0.2) # generate an bammdata object 
head(ed$eventData)


# Overall best shift configuration ----------------------------------------
best <- getBestShiftConfiguration(ed, expectedNumberOfShifts=1, threshold=5) # show the most frequent shift configuration in your credible shift set
plot.bammdata(best, lwd=1.2)
addBAMMshifts(best, cex=1.5)


# Extracting tip-specific evolutionary rates -----------------------------------------
tip.rates <- getTipRates(ed)
str(tip.rates)

speciation_mean<-tip.rates$lambda.avg # mean speciation rates
extinction_mean<-tip.rates$mu.avg # mean extinction rates
net_div<- speciation_mean-extinction_mean # net diversification rates (speciation minus extinction)

TipEstmates<-cbind.data.frame(Speciation=speciation_mean,Extinction=extinction_mean,NetDiv=net_div) # combining extracted data
head(TipEstmates)


# Exporting tip-specific evolutionary rates -------------------------------
# write.csv(TipEstmates,"BAMM_estimates_2638spp.csv")


# Histogram of the extracted evolutionary rates --------------------------
hist(TipEstmates$Speciation,xlab="average speciation") # mean speciation 
hist(TipEstmates$Extinction,xlab="average extinction") # mean extinction
hist(TipEstmates$NetDiv,xlab="average net diversification rates") # mean extinction


# Look at the speciation rates inferred -----------------------------------
bamm.fish <- plot.bammdata(ed, lwd=1, labels = F, spex = "s", method="polar",pal="Spectral")
addBAMMlegend(bamm.fish, nTicks = 4, side = 4, las = 1)


# Look at the extinction rates inferred -----------------------------------
bamm.fish <- plot.bammdata(ed, lwd=1, labels = F, spex = "e", method="polar",pal="Spectral")
addBAMMlegend(bamm.fish, nTicks = 4, side = 4, las = 1)


