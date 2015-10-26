# load libraries
# SESYNC and SLURM
setwd("stopover")
# devtools::install_github("SESYNC-ci/rslurm")

library(rslurm)
library(devtools)
library(msm)
library(parallel)

# install.packages("StopoverCode", repos = NULL, type="source")
library(StopoverCode)
# devtools::load_all("StopoverCode")
# this was replaced by 'StopoverCode' package, loaded with library
# source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

# load data
count_array <- readRDS('count_array_expanded.rds')
cov_array <- readRDS('covariates_array_expanded.rds')
cov_sites <- readRDS('covariates_sites_expanded.rds')
species_list <- readRDS('species_expanded.rds')

# choose parameter ranges
nBoots <- 100 # bootstraps to perform for each parameter combination
species <- 13 # corresponds to row in species_list
raw_cutoff <- c(5, 10, 20)
# Covariates for p (detection probability)
# 7 matrices, already scaled
# In order, temperature, temperature^2, wind, %cloud, survey duration, hour of day, #species recorded
p_cov1 <- 5 # Select detection covariates here (1:7 possible)
p_cov2 <- c(NA, 1:4, 6, 7) # Select detection covariates here (1:7 possible)
site_covs <- c("gdd", "lat") # for mu, w 
M <- c(1, 2, 3) #number of broods to estimate

params <- expand.grid(species, raw_cutoff, p_cov1, p_cov2, site_covs, M, stringsAsFactors = FALSE)
names(params) <- c("species", "raw_cutoff", "p_cov1", "p_cov2", "site_covs", "M")
params$nBoots <- nBoots

# data_file Rdata
dataIN <- c("count_array", "cov_array", "cov_sites", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))


# slurm.apply function
# now includes data prep, bootstrap site selection within function
# Question for Eleni: should same bootstrap sites be used for 
# different permutations of covariates to test them against each other?

SlurmStopover <- function(nRun){
  # parameters
  # is <<- necessary so inner functions can find all these?
  pars <- params[nRun, ]
  species <- pars$species
  raw_cutoff <- pars$raw_cutoff
  p_cov1 <- pars$p_cov1
  p_cov2 <- pars$p_cov2
  site_covs <- pars$site_covs
  M <<- pars$M
  nBoots <- pars$nBoots
  
  # these won't change, choose whether stopover model 
  # parameters have covariates or not
  p.m <<- "cov"
  w.m <<- "cov"
  mu.m <<- "cov"
  sigma.m <<- "het"
  phi.m <<- "const"
  
  ### data prep ###
  counts <- count_array[, , species] 
  
  # select sites with enough individuals counted
  siteRows <- which(rowSums(counts, na.rm = TRUE) >= raw_cutoff)
  counts <- counts[siteRows, ]
  counts[is.na(counts)] <- -1 #reassign NA for C program
  ALLcounts <- counts # ALL distinguishes from bootstrapped counts
  
  S <<- dim(counts)[1] 
  K <<- dim(counts)[2]  
  TIME <<- dim(counts)[2]
  
  # covariates
  # detection probability
  covs <- c(p_cov1, p_cov2)
  covs <- covs[!is.na(covs)]
  if (length(covs) > 1) cov.p <- cov_array[, , covs]
  if (length(covs) == 1) cov.p <- array(data = cov_array[, , covs], dim = c(dim(cov_array[,,covs]), 1))
  cov.p <- cov.p[siteRows, , , drop = FALSE]
  # qp is the number of covariates for p
  qp <- dim(cov.p)[3]
  for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
  cov.p[is.na(cov.p)] <- -1
  ALLcov.p <- cov.p
  
  # Time covariate for phi
  # Probably not necessary, since I say phi is constant
  cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))), S, K-1, byrow=TRUE) 
  ALLcov.phi <- cov.phi
  
  # Site covariates for mu and weight of broods
  if (site_covs == "gdd") s_cov <- scale(cov_sites$YearGDD)
  if (site_covs == "lat") s_cov <- scale(cov_sites$lat)
  s_cov <- s_cov[siteRows]
  ALLcov.w <- s_cov
  ALLcov.mu <- s_cov
  
  bootFits = list()
  # Bootstrap sites and fit model
  bootIdxs = array(0, c(nBoots, dim(counts)[1]))
  for (b in 1:nBoots){
    bSites <- sample(dim(counts)[1], replace = TRUE)
    startTime <- Sys.time()
    
    counts <<- ALLcounts[bSites, ]
    cov.w <<- ALLcov.w[bSites]
    cov.mu <<- ALLcov.mu[bSites]
    cov.p <<- ALLcov.p[bSites, , , drop = FALSE]
    qp <<- dim(cov.p)[3]
    
    ####
    Tries <- 3
    temp.fit <- list()
    temp.ll <- rep(NA, Tries)
    
    for (k in 1:Tries){
      start.list <<- startVals(p.m,w.m,mu.m,sigma.m,phi.m)
      pars.start <<- c(start.list$N,start.list$cvec, start.list$d0, start.list$d1, start.list$b0, start.list$b1,start.list$sigma,  start.list$a0, start.list$a1, start.list$a2) #this line remains the same for all models
      temp.fit[[k]] <- mLLMixtCounts.fit(p.m,w.m,mu.m,sigma.m,phi.m)
      temp.ll[k] <- temp.fit[[k]]$ll.val
    }
    
    if (length(which(is.na(temp.ll) == TRUE)) < Tries){
      tempchoose <- min(c(1:Tries)[which(temp.ll==max(temp.ll,na.rm=TRUE))])
      temp <-temp.fit[[tempchoose]]
    }else{
      temp <- list()
      temp$ll.val <- NA
    }
    
    temp$time <- startTime - Sys.time()
    temp$bSites <- bSites
    temp$pars <- pars
    temp$bootstrap <- b
    bootFits[[b]] <- temp
  }
  
  return(bootFits)
}

# test run of comparing different covariates with 3 bootstraps for each combo
sjob1 <- slurm_apply(f = SlurmStopover, params = paramIN, 
                     cpus_per_node = 8, nodes = 3, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

# what to do with output?
# first, want to see how many bootstraps ended in NA ll.val
# length of output list is nrow(param)*nBoots
outList <- get_slurm_out(sjob1)
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]]$pars
  out$bootstrap <- outList[[i]]$bootstrap
  out$ll.val <- outList[[i]]$ll.val
  out$time <- as.double(outList[[i]]$time, units = "mins")
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)

# test run of comparing different covariates with 3 bootstraps for each combo
sjob2 <- slurm_apply(f = SlurmStopover, params = paramIN, 
                     cpus_per_node = 8, nodes = 8, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

#5:39pm start 10/14

# what to do with output?
# first, want to see how many bootstraps ended in NA ll.val
# length of output list is nrow(param)*nBoots
outList <- get_slurm_out(sjob1)


