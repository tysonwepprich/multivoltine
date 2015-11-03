# Boostrap with a likelihood ratio test to determine number of generations automatically
# Steps:
# Fit model with M=1 (or other null hypothesis)
# Simulate bootstrap samples from MLE from model fit to null hypothesis
# Fit model to bootstrap samples for M=1 and M=2 (or other alternative hypothesis)
# Calculate log likelihood ratio statistic each time

# Computing approach
# Fit all potential null model(s) first using slurm cluster
# bootstrap simulations wrapped in function for cluster/parallel


# load libraries
# SESYNC and SLURM
# setwd("stopover")
# devtools::install_github("SESYNC-ci/rslurm")

# library(rslurm)
library(devtools)
library(msm)
library(parallel)
library(plyr)
library(dplyr)

# Eleni's functions
# this was replaced by 'StopoverCode' package, loaded with library
# source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

# for linux/sesync cluster
# remove.packages('StopoverCode') # do this to rebuild after edits
# install.packages("StopoverCode", repos = NULL, type="source")
# library(StopoverCode)

# for windows laptop
# load_all works, not install.package
devtools::load_all("StopoverCode", recompile = TRUE)

# load data
count_array <- readRDS('count_array_expanded.rds')
cov_array <- readRDS('covariates_array_expanded.rds')
cov_sites <- readRDS('covariates_sites_expanded.rds')
species_list <- readRDS('species_expanded.rds')
# 
# count_array <- readRDS('count_array.2009.rds')
# cov_array <- readRDS('covariates_array.2009.rds')
# cov_sites <- readRDS('covariates_sites.2009.rds')
# species_list_2009 <- readRDS('species.2009.rds')

# choose parameter ranges
species <- c(9:17) # corresponds to row in species_list
raw_cutoff <- 10 # c(5, 10, 20) higher cutoff increases fit, decreases comp. time
# Covariates for p (detection probability)
# 7 matrices, already scaled
# In order, temperature, temperature^2, wind, %cloud, survey duration, hour of day, #species recorded
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- "AnnGDD" # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(1:5) #number of broods to estimate

# can't just use these in expand.grid, because "common" limits options of other params
# what to do: find best covariate models first, then test against "common" 

# these 3 determined by "common" levels in p_cov and site_covs
# p.m <- c("cov", "common")[1]
# w.m <- c("cov", "common")[1]
# mu.m <- c("cov", "common")[1]
sigma.m <- "het" #  c("het", "hom")
# try phi now with survey "week", code up test of week vs. calendar day for age/time covariates
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")
# params$nBoots <- nBoots

# remove redundant parameter combos for detection, needed if p_cov1 and p_cov2 overlap
# params <- params[which(as.character(params$p_cov1) < as.character(params$p_cov2)), ]


# data_file Rdata
dataIN <- c("count_array", "cov_array", "cov_sites", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

# slurm.apply function
# fits each set of covariates once (multiple tries to find global LL)

SlurmCovs <- function(nRun){
  # parameters
  # is <<- necessary so inner functions can find all these?
  pars <- params[nRun, ]
  species <- pars$species
  raw_cutoff <- pars$raw_cutoff
  p_cov1 <- pars$p_cov1
  p_cov2 <- pars$p_cov2
  site_covs <- pars$site_covs
  M <<- pars$M
  #   p.m <<- pars$p.m
  #   w.m <<- pars$w.m
  #   mu.m <<- pars$mu.m
  sigma.m <<- pars$sigma.m
  phi.m <<- pars$phi.m
  
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
  covs <- sort(as.numeric(covs[covs %in% as.character(c(1:7))]))
  if (length(covs) == 0) {
    p.m <<- "common" 
    cov.p <- NA
  } else {
    p.m <<- "cov"
    if (length(covs) > 1) cov.p <- cov_array[, , covs]
    if (length(covs) == 1) cov.p <- array(data = cov_array[, , covs], dim = c(dim(cov_array[,,covs]), 1))
    cov.p <- cov.p[siteRows, , , drop = FALSE]
    # qp is the number of covariates for p
    qp <- dim(cov.p)[3]
    for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
    cov.p[is.na(cov.p)] <- -1
  }
  ALLcov.p <- cov.p
  
  
  # Time (week) covariate for phi
  cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))), S, K-1, byrow=TRUE) 
  ALLcov.phi <- cov.phi
  
  # Site covariates for mu and weight of broods
  if (site_covs == "common") {
    s_cov <- NA
    mu.m <<- "common"
    w.m <<- "common"
  }else{
    mu.m <<- "cov"
    w.m <<- "cov"
    if (site_covs == "AnnGDD") s_cov <- scale(cov_sites$YearGDD)
    if (site_covs == "SprGDD") s_cov <- scale(cov_sites$SpringGDD)
    if (site_covs == "lat") s_cov <- scale(cov_sites$lat)
    s_cov <- s_cov[siteRows]
  }
  ALLcov.w <- s_cov
  ALLcov.mu <- s_cov
  
  startTime <- Sys.time()
  
  counts <<- ALLcounts
  cov.w <<- ALLcov.w
  cov.mu <<- ALLcov.mu
  cov.p <<- ALLcov.p
  cov.phi <<- ALLcov.phi
  qp <<- dim(cov.p)[3]
  
  ####
  out <- list()
  Tries <- 5
  temp.fit <- list()
  temp.ll <- rep(NA, Tries)
  
  for (k in 1:Tries){
    start.list <<- startVals(p.m,w.m,mu.m,sigma.m,phi.m)
    pars.start <<- c(start.list$N,start.list$cvec, start.list$d0, start.list$d1, start.list$b0, start.list$b1,start.list$sigma,  start.list$a0, start.list$a1, start.list$a2) #this line remains the same for all models
    temp.fit[[k]] <- mLLMixtCounts.fit(p.m,w.m,mu.m,sigma.m,phi.m)
    temp.ll[k] <- temp.fit[[k]]$ll.val
  }
  
  if (length(which(is.na(temp.ll) == TRUE)) < Tries){
    tempchoose <- min(c(1:Tries)[which(temp.ll==max(temp.ll, na.rm=TRUE))])
    temp <-temp.fit[[tempchoose]]
  }else{
    temp <- list()
    temp$ll.val <- NA
  }
  
  temp$time <- startTime - Sys.time()
  temp$pars <- pars
  out[[1]] <- temp
  return(out)
}


# cl <- makeCluster(3)
# clusterEvalQ(cl, {
#   library(devtools)
#   library(msm)
#   library(dplyr)
#   devtools::load_all("StopoverCode", recompile = TRUE)
#   load("dataIN.RData")
# })
# test <- parLapply(cl, paramIN$nRun, SlurmCovs)
# stopCluster(cl)
# 



# calculate null hypotheses for M = 1:5 for different species
baseline <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 4, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")




