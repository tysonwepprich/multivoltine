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

library(rslurm)
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
library(StopoverCode)

# for windows laptop
# load_all works, not install.package
# devtools::load_all("StopoverCode", recompile = TRUE)

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
M <- c(1, 1, 1) # c(1:5) #number of broods to estimate

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
                     cpus_per_node = 8, nodes = 6, 
                     data_file = "dataIN.RData", 
                     # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                     output = "raw")

# calculate null hypotheses for M = 1 for different species
# got a lot of LL.val = NA for M = 1 the first time through, limits p-value
# this is just to test if more tries means its possible to fit
baseline1M <- slurm_apply(f = SlurmCovs, params = paramIN, 
                        cpus_per_node = 8, nodes = 2, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")

# simulation of sample from MLE parameters for each species/M hypothesis

# slurm_out <- get_slurm_out(baseline)
# saveRDS(slurm_out, "baselineMtest.rds")
base_out <- readRDS("baselineMtest.rds")
# base1_out <- get_slurm_out(baseline1M)

outList <- base_out
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
baselineDF <- outDF

# what to do about ll.val == NA for M = 1?
nsim <- 100
SampleList <- vector("list", length = length(slurm_out)*nsim)
# make bootstrap data simulations for each model x 100
for (bs in 1:length(SampleList)){
  mod <- (bs + nsim - 1) %/% nsim        #ADD IN INDEX FOR MOD AND BOOTSTRAP TO TRACK THESE
  nullFit <- slurm_out[[mod]]
  # building output list
  SampleList[[bs]]$pars <- nullFit$pars
  # move on to next model if ll.val is NA
  if (is.na(nullFit$ll.val)){
    SampleList[[bs]]$ll.val <- NA
    SampleList[[bs]]$npar <- NA
    next
  }else{
    SampleList[[bs]]$ll.val <- nullFit$ll.val
    SampleList[[bs]]$npar <- nullFit$npar
  }
  
  species <- nullFit$pars$species
  M <- nullFit$pars$M
  counts <- count_array[, , species] 
  
  # select sites with enough individuals counted
  siteRows <- which(rowSums(counts, na.rm = TRUE) >= raw_cutoff)
  counts <- counts[siteRows, ]
  
  S <<- dim(counts)[1] 
  K <<- dim(counts)[2]  
  TIME <<- dim(counts)[2]
  
  N.est <- nullFit$N.est
  phi.est <- nullFit$phi.est[1,1,1]
  c0.est <- nullFit$cest[1]
  c1.est <- nullFit$cest[2]
  b0.est <- nullFit$b0.est
  b1.est <- nullFit$b1.est
  d0.est <- nullFit$d0.est
  d1.est <- nullFit$d1.est
  sigma.est <- nullFit$sigma.est
  w.est <- nullFit$w.est
  
  simData <- list()
    N.tr <- rpois(S, N.est)
    phi.tr <- matrix(phi.est, TIME-1, TIME-1)
    
    cov.p_vary <- array(runif(S*TIME, -1, 1), c(S, TIME))
    
    c0.tr <- c0.est
    c1.tr <- c1.est
    p_vary.tr <- expo(c0.tr + c1.tr*cov.p_vary)
    
    cov.site_all <- rnorm(S)
    b0.tr <- matrix(rep(b0.est, S), ncol = M, byrow = TRUE)
    b1.tr <- b1.est
    mu.tr <- exp(b0.tr + b1.tr*cov.site_all)
    d0.tr <- d0.est
    d1.tr <- d1.est
    sd.tr <- sigma.est
    
    if (ncol(w.est) == 1){
      w.tr <- w.est
    }else{
      w.tr <- matrix(NA, ncol = M, nrow = S)
      for(site in 1:S){
        w.tr[site,] <- lgp(d0.tr + d1.tr * cov.site_all[site]) 
      }    
    }
    
    betta.tr <- matrix(0, S, TIME)
    if (M == 1){
      for(s in 1:S){
        betta.tr[s,] <-  c(pnorm(1, mean=mu.tr[s], sd=sd.tr[s]),
                           pnorm(2:(TIME-1),mean=mu.tr[s],sd=sd.tr[s])-pnorm(1:(TIME-2), mean=mu.tr[s], sd=sd.tr[s]),
                           1-pnorm(TIME-1,mean=mu.tr[s],sd=sd.tr[s]))
      }
    }else{
      for(s in 1:S){
        betta.site <- rep(0,TIME)
        for(m in 1:M){
          betta.site <- betta.site + c(w.tr[s,m]*c(pnorm(1,mean=mu.tr[s,m],sd=sd.tr[s,m]),
                                                   pnorm(2:(TIME-1), mean=mu.tr[s,m], sd=sd.tr[s,m])-pnorm(1:(TIME-2),
                                                                                                           mean=mu.tr[s,m],sd=sd.tr[s,m]), 1-pnorm(TIME-1,mean=mu.tr[s,m],sd=sd.tr[s,m])))
        }
        betta.tr[s,] <- betta.site
      }
    }
    
    
    
    lambda.tr <- array(0, c(S, TIME))
    counts_all <- array(0,c(S, TIME))
    
    for(i in 1:S){
      
      lambda.tr[i,] <- N.tr[i]*betta.tr[i,]
      counts_all[i,] <- rpois(TIME, lambda.tr[i,]*p_vary.tr[i,])	
      
      for(j in 2:TIME){
        for(b in 1:(j-1)){
          lambda.tr[i,j] <- lambda.tr[i,j] + N.tr[i]*betta.tr[i,b]*prod(phi.tr[b,b:(j-1)])
        }
        counts_all[i,j] <- rpois(1,lambda.tr[i,j]*p_vary.tr[i,j])
      }
    }
    #drop records to approximate rate of missing survey each week
    
    counts_missing <- array(0, c(S, TIME))
    
    PropSurv <- function(vec){
      prop <- length(which(is.na(vec) == FALSE))/length(vec)
      return(prop)
    }
    FreqSurv <- apply(counts, 2, PropSurv)
    for (s in 1:S){
      counts_missing[s,] <- rbinom(TIME, 1, FreqSurv)
    }
    counts_missing[counts_missing == 0] <- NA
    counts_missing <- counts_all * counts_missing
    counts_missing[is.na(counts_missing)] <- -1 #reassign NA for C program
    
    
    simData$counts <- counts_missing
    simData$p_cov <- cov.p_vary
    simData$site_cov <- cov.site_all

    SampleList[[bs]]$simData <- simData
} #close Sample for loop

# saveRDS(SampleList, file = "SampleList.rds")
SampleList <- readRDS("SampleList.rds")

# fit the bootstrap sample datasets using SLURM

# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:length(SampleList)))

# Fit models for null and alternative (M + 1) for each simulated dataset
SlurmGeneration <- function(nRun){
  out <- vector("list", length = 2)
  BSdata <- SampleList[[nRun]]
  if (is.na(BSdata$ll.val)){
    out[[1]]$pars <- BSdata$pars
    out[[1]]$ll.val <- NA
    out[[1]]$nRun <- nRun
    out[[2]] <- NULL
  }else{
    # parameters
    pars <- BSdata$pars
    species <- pars$species
    raw_cutoff <- pars$raw_cutoff
    p_cov1 <- pars$p_cov1
    p_cov2 <- pars$p_cov2
    site_covs <- pars$site_covs
    M <<- pars$M
    sigma.m <<- pars$sigma.m
    phi.m <<- pars$phi.m
    p.m <<- "cov"
    w.m <<- "cov"
    mu.m <<- "cov"
    ### data prep ###
    counts <<- BSdata$simData$counts
    
    S <<- dim(counts)[1] 
    K <<- dim(counts)[2]  
    TIME <<- dim(counts)[2]
    
    p_cov <- BSdata$simData$p_cov
    p_cov[which(counts == -1)] <- NA
    
    cov.p <- array(data = p_cov, dim = c(dim(p_cov), 1))
    qp <- dim(cov.p)[3]
    for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
    cov.p[is.na(cov.p)] <- -1
    ALLcov.p <- cov.p
    
    # Time (week) covariate for phi
    cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))), S, K-1, byrow=TRUE) 
    ALLcov.phi <- cov.phi
    
    s_cov <- scale(BSdata$simData$site_cov)
    ALLcov.w <- s_cov
    ALLcov.mu <- s_cov
    
    cov.w <<- ALLcov.w
    cov.mu <<- ALLcov.mu
    cov.p <<- ALLcov.p
    cov.phi <<- ALLcov.phi
    qp <<- dim(cov.p)[3]
    
    
    ####
    # Fit null hypothesis for M
    Tries <- 5
    temp.fit <- list()
    temp.ll <- rep(NA, Tries)
    startTime <- Sys.time()
    
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
    temp$model <- "null"
    temp$time <- startTime - Sys.time()
    temp$pars <- pars
    out[[1]] <- temp
    out[[1]]$nRun <- nRun
    
    # Fit alternative model
    M <<- M + 1
    Tries <- 5
    temp.fit <- list()
    temp.ll <- rep(NA, Tries)
    startTime <- Sys.time()
    
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
    temp$model <- "alt"
    temp$time <- startTime - Sys.time()
    temp$pars <- pars
    out[[2]] <- temp
    out[[2]]$nRun <- nRun
  } #close ifelse
  return(out)
}



# calculate null hypotheses for M = 1:5 for different species
alt <- slurm_apply(f = SlurmGeneration, params = paramIN, 
                        cpus_per_node = 8, nodes = 12, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")


test2 <- get_slurm_out(alt)
# finding results when you clear the workspace and get_slurm_out no longer works
slurm_codes <- c("slr2747")
slurm_out <- list()

for (j in 1:length(slurm_codes)){
  missing_files <- c()
  tmpEnv <- new.env()
  for (i in 0:11) {
    fname <- paste0(slurm_codes[j], "_", i, 
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
}

test <- do.call(rbind, lapply(slurm_out, function(x) length(x)))
####PROBLEM: output list has 2 components, each needed to be collapsed into a row in outDF
# many cases where M=1 doesn't work but M=2 does for same bootstrap sample

####NEW PROBLEM: some output in list of 1, but model still fit for null and alt, but retrieved
# in 2 separate lists. Now have slurm_out elements with length 1, 2, 5 or 22!

outList <- slurm_out
outDF <- data.frame()
for (i in 1:length(outList)){
  # expected format length 2
  if (length(outList[[i]]) %in% c(1,2)){
    for (j in 1:length(outList[[i]])){
      species <- outList[[i]][[j]]$pars$species
      nRun <- outList[[i]][[j]]$nRun
      ll.val <- outList[[i]][[j]]$ll.val
      if (is.na(ll.val)){
        npar <- NA
        M <- outList[[i]][[j]]$pars$M
        model <- NA
      }else{
        npar <- outList[[i]][[j]]$npar
        M <- dim(outList[[i]][[j]]$mu.est)[2]
        model <- outList[[i]][[j]]$model
      }
      outDF <- rbind(outDF, data.frame(species, nRun, M, model, ll.val, npar))
    }
  }
  if (length(outList[[i]]) %in% c(5, 22)){
    species <- outList[[i]]$pars$species
    nRun <- outList[[i]]$nRun
    ll.val <- outList[[i]]$ll.val
    if (is.na(ll.val)){
      npar <- NA
      M <- outList[[i]]$pars$M
      model <- NA
    }else{
      npar <- outList[[i]]$npar
      M <- dim(outList[[i]]$mu.est)[2]
      model <- outList[[i]]$model
    }
    outDF <- rbind(outDF, data.frame(species, nRun, M, model, ll.val, npar))
  }
}
BSmods <- outDF



# from output, calculate 
# -2loglambda for each bootstrap data sample
# compare to -2loglambda from models fit to all data
# to derive p-value based on
# alpha = 1 - j / (B+1)


BSpval <- function(nullM, spec){
  nullM <- nullM
  spec <- spec
  origNull <- baselineDF %>% filter(species == spec, M == nullM) %>% select(ll.val)
  origAlt <- baselineDF %>% filter(species == spec, M == nullM + 1) %>% select(ll.val)
  # problem when null mod generally doesn't fit, many alternative mods do
  
  BStests <- BSmods %>% filter(species == spec) %>% 
    filter(M == nullM & model == "null" | M == nullM + 1 & model == "alt") %>%
    group_by(nRun) %>% mutate(NumSuccess = length(ll.val)) %>%
    filter(NumSuccess == 2)
  
#   if (length(which(BStests$model == "null")) < length(which(BStests$model == "alt"))){
#     
#   }
  
  LRdistr <- BStests %>% group_by(nRun) %>%
    summarise(neg2logLam = 2 * (ll.val[model == "alt"] - ll.val[model == "null"]))
  B <- nrow(LRdistr)
  # alpha <- 0.05
  origLR <- 2 * (origAlt - origNull)
  j <- length(which(LRdistr$neg2logLam < origLR$ll.val)) 
  # P <- 1 - (3 * j - 1) / (3 * B + 1)
  P <- 1 - j / (B + 1)
  return(P)
}
# if P > alpha, alternative hypothesis is better

# species 11, getting some runs where null has slightly lower deviance than alternative
# this causes negative values of neg2logLam, which should be impossible

expand.grid.alt <- function(seq1,seq2) {
  cbind(rep.int(seq1, length(seq2)),
        c(t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
}

parsIN <- data.frame(expand.grid.alt(c(1:4), unique(BSmods$species)))
names(parsIN) <- c("nullM", "spec")
Mtest <- parsIN %>% rowwise() %>% mutate(pval = BSpval(nullM, spec)) %>% data.frame()
# saveRDS(Mtest, file = "MtestTry1.rds")

