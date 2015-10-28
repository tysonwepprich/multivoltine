# Update slurmRum.R to better parallelize bootstrap
# 1st try puts all bootstraps on one CPU, each CPU a unique BS run of param combo
# Now, make slurm.apply function run a single bootstrap, should be more efficient

# Also, this file will not run bootstrap for different covariates combos
# Use AIC like in Eleni's paper to assess M and covariates
# Add covariates for phi like age and time, since they were important in her paper.

# Run covariate test for each species, multiple years?

# load libraries
# SESYNC and SLURM
# setwd("stopover")
# devtools::install_github("SESYNC-ci/rslurm")

library(rslurm)
library(devtools)
library(msm)
library(parallel)
library(dplyr)
# on Windows laptop, load_all works, not install.package
# install.packages("StopoverCode", repos = NULL, type="source")
# library(StopoverCode)
devtools::load_all("StopoverCode", recompile = FALSE) # recompile true when first time loading on windows
# this was replaced by 'StopoverCode' package, loaded with library
# source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

# load data
count_array <- readRDS('count_array_expanded.rds')
cov_array <- readRDS('covariates_array_expanded.rds')
cov_sites <- readRDS('covariates_sites_expanded.rds')
species_list <- readRDS('species_expanded.rds')

# choose parameter ranges
# nBoots <- 100 # bootstraps to perform for each parameter combination
species <- c(10, 11, 12) # corresponds to row in species_list
raw_cutoff <- 10 # c(5, 10, 20) higher cutoff increases fit, decreases comp. time
# Covariates for p (detection probability)
# 7 matrices, already scaled
# In order, temperature, temperature^2, wind, %cloud, survey duration, hour of day, #species recorded
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(1, 2, 3) #number of broods to estimate
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- "AnnGDD" #c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(2) #number of broods to estimate

# can't just use these in expand.grid, because "common" limits options of other params
# what to do: find best covariate models first, then test against "common" 

# these 3 determined by "common" levels in p_cov and site_covs
# p.m <- c("cov", "common")[1]
# w.m <- c("cov", "common")[1]
# mu.m <- c("cov", "common")[1]
sigma.m <- "het" # c("het", "hom")
# try phi now with survey "week", code up test of week vs. calendar day for age/time covariates
phi.m <- c("const", "logit.a")

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
# now includes data prep, bootstrap site selection within function
# Question for Eleni: should same bootstrap sites be used for 
# different permutations of covariates to test them against each other?

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

cl <- makeCluster(3)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  devtools::load_all("StopoverCode", recompile = TRUE)
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmCovs)
stopCluster(cl)





# test run of comparing different covariates with 3 bootstraps for each combo
sjob1 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 10, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

# what to do with output?
# first, want to see how many bootstraps ended in NA ll.val
# length of output list is nrow(param)*nBoots
jobs <- list(sjob1, sjob2, sjob3, sjob4)
jobs <- list(sjob5, sjob6, sjob7)
jobs <- list(sjob8)
outALL <- data.frame()
for (j in 1:length(jobs)){
  outList <- get_slurm_out(jobs[[j]])
  outDF <- list()
  for (i in 1:length(outList)){
    out <- outList[[i]]$pars
    out$model <- i
    out$ll.val <- outList[[i]]$ll.val
    if (is.na(out$ll.val)){
      out$npar <- NA
    }else{
      out$npar <- outList[[i]]$npar
    }
    out$time <- as.double(outList[[i]]$time, units = "mins")
    outDF[[i]] <- out
  }
  outDF <- do.call("rbind", outDF)
  outALL <- rbind(outALL, outDF)
}
  
# something wrong with phi covariates, only gives NA ll.val if not 'const'

covTest <- outALL[is.na(outALL$ll.val) == FALSE, ]
covTest$AIC <- -2 * covTest$ll.val + 2 * covTest$npar
# covTest$AICc <- -2 * covTest$ll.val + 2 * covTest$npar
# minAIC <- min(covTest$AIC)
# covTest$weight1 <- exp(-0.5 * (covTest$AIC - minAIC))
# sumweight <- sum(covTest$weight1)
# covTest$weight <- covTest$weight1 / sumweight
  
test <- covTest %>%
  group_by(species) %>%
  mutate(weight1 = exp(-0.5 * (AIC - min(AIC)))) %>%
  mutate(weight = weight1 / sum(weight1))

test2 <- test %>%
  select(-raw_cutoff, -weight1) %>%
  filter(species == 10) %>%
  arrange(AIC) %>%
  data.frame()

# tested 4 species with covariates, more M than expected...
# list-length detection covariates best for some, maybe duration isn't needed?

# extract models with different phi structures to visualize
phitest <- test2[test2$site_covs == "lat" & test2$p_cov2 == 3 & test2$sigma.m == "het",]
phimods <- phitest$model[1:5]

phiout <- get_slurm_out(sjob8)
philist <- phiout[phimods]
saveRDS(philist, "philist.RDS")

# new function edits to allow input params to have "common" covariate specifications
# paramIN <- data.frame(nRun = seq(1:773))
sjob5 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 4, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

# paramIN <- data.frame(nRun = seq(774:1546))
sjob6 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 4, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")


# paramIN <- data.frame(nRun = seq(1547:nrow(params)))
sjob7 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 4, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")


# rerun, holding p_cov1 as list-length
# species 10
sjob8 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 8, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

# species 11
sjob9 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 8, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")


# species 12
sjob10 <- slurm_apply(f = SlurmCovs, params = paramIN, 
                     cpus_per_node = 8, nodes = 8, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")


# finding results when you clear the workspace and get_slurm_out no longer works

slurm_out <- list()
missing_files <- c()
tmpEnv <- new.env()
for (i in 0:7) {
  fname <- paste0("slr8096", "_", i, 
                  ".RData")
  if (fname %in% dir()) {
    load(fname, envir = tmpEnv)
    slurm_out <- c(slurm_out, get(".rslurm_result", 
                                  envir = tmpEnv))
  }
  else {
    missing_files <- c(missing_files, fname)
  }
}

outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]]$pars
  out$model <- i
  out$ll.val <- outList[[i]]$ll.val
  if (is.na(out$ll.val)){
    out$npar <- NA
  }else{
    out$npar <- outList[[i]]$npar
  }
  out$time <- as.double(outList[[i]]$time, units = "mins")
  outDF[[i]] <- out
}
outALL <- do.call("rbind", outDF)

