# source this for bootstrap M functions

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


# expand.grid altered for integers, not just factors/characters
expand.grid.alt <- function(seq1,seq2) {
  cbind(rep.int(seq1, length(seq2)),
        c(t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
}
