source('bootstrapMfunctions.R')

allSpecies <- read.csv("data/MultivoltineSpecies.csv", header = TRUE)
i <- 4
species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i] + 1

dat <- SpeciesDataP1(species)

library(PReMiuM)

inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
                     xModel=inputs$xModel, nSweeps=10, nClusInit=20,
                     nBurn=20, data=inputs$inputData, output="output", 
                     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                     fixedEffectsNames = inputs$fixedEffectNames)


library(mclust)
data(acidity)
mod4 = densityMclust(acidity)
summary(mod4)



library(flexmix)
# this one seems promising, allows covariates, repeated measures.