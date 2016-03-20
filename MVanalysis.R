# Hopefully final draft of stopover model analysis 
# 1. By species, constrain M to best values (by manually looking at phenograms)
# 2. Auto-bootstrap M selection for base M and M+1 for each year
# 3. Run bootstrap to get parameter estimates for best M for each year

# Difficult thing is optimizing for slurm. May not be able to have one big function,
# since some steps may need many cores, some few

# file needs to be edited depending on whether Windows or Linux for package loading
source('bootstrapMfunctions.R')

allSpecies <- read.csv("data/MultivoltineSpecies.csv", header = TRUE)

# choose your species
# running all species with minimum number of covariates to test M
# hoped that this would make absurdly high populations less likely
# even some still have them, maybe run them over again and again until they find a better fit
# threshold: either > 1000, or maybe if same year, different M give wildly varying N estimates (5x or 10x different)


# 14 RSP 6698 ** 9718 ** 2241 GEN
# 18 WIDW 3540 ** 2377
# 16 Spice 3941
# 19 Zab 1497 ** 5423
# 17 Vice 1621
# 13 Peck 1727 ** 5312
# 1 Black Swal 1933 ** 4162
# 2 CWN 2020 **4282
# 3 ETS 2153 ** 1668
# 4 Euro 2285
# 5 Hack 2493
# 6 Hobo 2599
# 7 Juv 2864
# 8 Least 2967
# 15 SSSkip 3152
# 9 LGW 3308
# 10 LWS 3441 ** 476 ** 5332 GEN
# 11 NBD 3528
# 12 NPE 3629
i <- 4

species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i] + 1

# somewhat unwieldly, list with each year as a list of 4 (year, counts, surv_covs, site_covs)
# dat <- SpeciesData(species)
dat <- SpeciesDataP1(species)


# for each species, select parameters
# how much data available for modeling?

count_cutoff <- 5
surv_cutoff <- 3
data_avail <- data.frame()
for (j in 1:length(dat)){
  list_index <- j
  temp <- dat[[j]]
  year <- temp[[1]]
  counts <- temp[[2]]
  sum_count <- apply(counts, 1, sum, na.rm = TRUE)
  surv_present <- apply(counts, 1, function(x) length(which(x > 0)))
  cutoff_met <- length(which(sum_count >= count_cutoff))
  survs_met <- length(which(surv_present >= surv_cutoff))
  both <- cbind(sum_count, surv_present)
  both_met <- length(which(apply(both, 1, function(x, count_cutoff, surv_cutoff) x[1] >= count_cutoff & 
                                   x[2] >= surv_cutoff, count_cutoff = count_cutoff, surv_cutoff = surv_cutoff)))
  new_row <- data.frame(count_cutoff, surv_cutoff, year, list_index, cutoff_met, survs_met, both_met)
  data_avail <- rbind(data_avail, new_row)
}


list_index_min_data <- unique(data_avail$list_index[data_avail$both_met >= 5])
# list_index_min_data <- c(6, 9) # test problem with SSSkip 


# choose parameter ranges
raw_cutoff <- 5 # c(5, 10)
p_cov1 <- "none" #c("none", 7) # Select detection covariates here (1:7 possible)
p_cov2 <- "none" #c("none", 1) # Select detection covariates here (1:7 possible)
p_cov3 <- "none" #c("none", 2)
p_cov4 <- "none" #c("none", 1)
site_covs <- "AnnGDD" # c("AnnGDD", "lat") # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(minBrood:maxBrood) #number of broods to estimate
sigma.m <- "hom" #c("het", "hom")
phi.m <- "const" #c("const", "quad.t")

params <- expand.grid(species, list_index_min_data, raw_cutoff, p_cov1, p_cov2, p_cov3, p_cov4,
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "list_index", "raw_cutoff", "p_cov1", "p_cov2", "p_cov3", "p_cov4",
                   "site_covs", "M", "sigma.m", "phi.m")
# double size of params with temperature p covariates
# params2 <- params
# params2$p_cov2 <- 1
# params2$p_cov3 <- 2
# params <- rbind(params, params2)
params <- rbind(params, params, params, params)

params$param_row <- 1:nrow(params)
params <- params[sample(1:nrow(params)), ] #rearrange for parallel comp speed

# data_file Rdata
dataIN5 <- c("dat", "params")
save(list = dataIN5, file = "dataIN5.RData")

# simple param file for slurm.apply
paramIN5 <- data.frame(nRun = seq(1:nrow(params)))

# # single core
# system.time({
#   test <- lapply(paramIN$nRun, SlurmCovs)
#   })


# multicore
system.time({
cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  library(StopoverCode) #on linux
  # devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapplyLB(cl, paramIN$nRun, SlurmCovs)
stopCluster(cl)
})

saveRDS(test, file = "EUROcovs.rds")

# calculate null hypotheses for M for different species
testCovs <- slurm_apply(f = SlurmCovs, params = paramIN, 
                          cpus_per_node = 8, nodes = 4, 
                          data_file = "dataIN.RData", 
                          # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                          output = "raw")


# Next step, go through processSlurmCov.R to choose best site_cov for species over all years
# May need to edit code depending on how lists are output, lengths may differ
# this throws indexing off for some species/years
# Stupidly manual, but at least allows for some checks.

# for a given species
# simulate data for each year
# fit M and M+1 to compare and get p-value to find # of mixture model modes

# I think calculating the MLE of MaxBrood + 1 was unnecssary for the following bootstrap
# could still use it to compare mode selection with AIC vs bootstrap though

# 100 simulations sufficient to see if p-value close to cutoff 0.05 or 0.1
# might run more if near boundary

# need to get slurm_out2 for each species before simulation

########################
# extract data from SlurmCov results from parlapply (non-slurm)

# results <- list.files("slurmCovOutput/otherResults/")
# for (res in 2:14){
# setwd("slurmCovOutput/otherResults/")
# temp <- readRDS(results[res])
temp <- readRDS("leastskipCOV.rds")
test <- do.call(rbind, lapply(temp, function(x) length(x[[1]]))) # extra layer of list


outList <- temp
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]][[1]]$pars
  out$model <- i
  out$ll.val <- outList[[i]][[1]]$ll.val
  if (is.na(out$ll.val)){
    out$npar <- NA
    out$maxNest <- NA
    out$medP <- NA
  }else{
    out$npar <- outList[[i]][[1]]$npar
    out$maxNest <- round(max(outList[[i]][[1]]$N.est))
    out$medP <- round(median(outList[[i]][[1]]$p.est), 3)
  }
  out$time <- as.double(outList[[i]][[1]]$time, units = "mins")
  out$nRun <- outList[[i]][[1]]$nRun
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
baselineDF <- outDF
# setwd("../../")

#############################################


# extract data from SlurmCov results
slurm_codes <- c("slr2888")
slurm_out <- list()
# setwd("slurmCovOutput")

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
# setwd("../")

saveRDS(slurm_out, "CWNslurmcovs.rds")


# slurm_out<- readRDS("RSPslurmcovs.rds")
outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  if (length(outList[[i]]) == 1){
    out <- NA
  }else{
    out <- outList[[i]]$pars
    out$model <- i
    out$ll.val <- outList[[i]]$ll.val
    if (is.na(out$ll.val)){
      out$npar <- NA
      out$maxNest <- NA
      out$medP <- NA

    }else{
      out$npar <- outList[[i]]$npar
      out$maxNest <- round(max(outList[[i]]$N.est))
      out$medP <- round(median(outList[[i]]$p.est), 3)
    }
    out$time <- as.double(outList[[i]]$time, units = "mins")
    out$nRun <- outList[[i]]$nRun
  }
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
baselineDF <- outDF

##################################################
# when run multiple time, variation in loglik and parameters for same M
# choose best model, then simulate data and compare to M+1 for LRT to choose M

# bestmods will be used for LRT statistic to compare to null distribution from simulations
bestmods <- baselineDF %>% group_by(list_index, M) %>%
  filter(ll.val == max(ll.val, na.rm = TRUE)) %>%
  filter(maxNest < 1000 * median(baselineDF$maxNest)) %>%    ##trying to deal with huge N estimates
  arrange(list_index)

# nullmods are model fits from which data is simulated 
nullmods <- bestmods %>% group_by(list_index) %>% filter(M < max(M))

# simFits <- slurm_out[c(nullmods$model)]
simFits <- temp[c(nullmods$model)]

###############################################
# get simulated data from null hypotheses for M (# of generations in a year)

# use only 50 nsims to start
SampleList <- SimNullData(simFits, nsim = 50, spec_data = dat)


saveRDS(SampleList, file = paste("simDataGenMode/", gsub(" ", "", species, fixed = TRUE), "simdata.rds", sep = ""))


# SampleList <- readRDS("simDataGenMode/HobomokSkippersimdata.rds")
# data_file Rdata

dataIN1 <- c("SampleList")
save(list = dataIN1, file = "dataIN1.RData")

# simple param file for slurm.apply
paramIN1 <- data.frame(nRun = sample(seq(1:length(SampleList)))) # random nRun so split even for parallel


cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  library(StopoverCode) #on linux
  # devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
time <- system.time({test <- parLapply(cl, paramIN$nRun, SlurmGenerationP1)})
stopCluster(cl)

saveRDS(test, file = "eurpBSmod.rds")

# calculate null hypotheses for same species, different years
peckBS <- slurm_apply(f = SlurmGenerationP1, params = paramIN1, 
                   cpus_per_node = 8, nodes = 3, 
                   data_file = "dataIN1.RData", 
                   # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                   output = "raw")





# extract data from SlurmGeneration results
slurm_codes <- c("slr5332")
slurm_out <- list()
# setwd("slurmCovOutput")

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
# setwd("../")

# outList <- readRDS("npeBSmod.rds")
outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  if (length(outList[[i]]) == 1){
    out <- NA
  }else{
    for (j in 1:length(outList[[i]])){
      out <- outList[[i]][[j]]$pars
      out$model <- outList[[i]][[j]]$model
      out$ll.val <- outList[[i]][[j]]$ll.val
      if (is.na(out$ll.val)){
        out$npar <- NA
        out$maxNest <- NA
        out$medP <- NA
        out$obsM <- NA
        
      }else{
        out$npar <- outList[[i]][[j]]$npar
        out$maxNest <- round(max(outList[[i]][[j]]$N.est))
        out$medP <- round(median(outList[[i]][[j]]$p.est), 3)
        out$obsM <- dim(outList[[i]][[j]]$mu.est)[2]
      }
      out$time <- as.double(outList[[i]][[j]]$time, units = "mins")
      out$nRun <- outList[[i]][[j]]$nRun
      outDF[[(length(outDF) + 1)]] <- out
    }
  }
}

outDF <- do.call("rbind", outDF)
BSmods <- outDF



# M's not accurate
BSmods$M[BSmods$model == "alt"] <- BSmods$M[BSmods$model == "alt"] + 1
# BSmods <- BSmods %>% filter(param_row == 12)

# from baselineDF from slurmCov
slurm_out<- readRDS("RSPslurmcovs.rds")

outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  if (length(outList[[i]]) == 1){
    out <- NA
  }else{
    out <- outList[[i]]$pars
    out$model <- i
    out$ll.val <- outList[[i]]$ll.val
    if (is.na(out$ll.val)){
      out$npar <- NA
      out$maxNest <- NA
      out$medP <- NA
      
    }else{
      out$npar <- outList[[i]]$npar
      out$maxNest <- round(max(outList[[i]]$N.est))
      out$medP <- round(median(outList[[i]]$p.est), 3)
    }
    out$time <- as.double(outList[[i]]$time, units = "mins")
    out$nRun <- outList[[i]]$nRun
  }
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
baselineDF <- outDF


bestmods <- baselineDF %>% group_by(list_index, M) %>%
  filter(ll.val == max(ll.val, na.rm = TRUE)) %>%
  arrange(list_index)


# tests <- baselineDF[c(2, 3, 7),]
# baselineDF <- tests[c(3, 1), ]
# baselineDF <- tests[c(3, 2), ]

listP <- list()
for (i in 1:length(unique(bestmods$list_index))){
  tempindex <- unique(bestmods$list_index)[i]
  tempbaselineDF <- bestmods %>% filter(list_index == tempindex)
  testM <- unique(tempbaselineDF$M)
  testM <- testM[-which(testM == max(testM))]
  out <- data.frame()
  for (m in testM){
    tempBSmods <- BSmods %>% filter(list_index == tempindex) %>%
      filter(M == m & model == "null" | M == m + 1 & model == "alt")
    if (nrow(tempBSmods) == 0){
      pval <- data.frame(nullM = m, pval = NA, nfits = 0, nullNA = NA, altNA = NA)
    } else {
      pval <- BSpval(nullM = m, spec = species, BSmods = tempBSmods, baselineDF = tempbaselineDF) # also needs BSmods from fitting simulated data, and baselineDF from original fitting
    }
    out <- rbind(out, pval)
  }
  out$list_index <- tempindex
  listP[[i]] <- out
}
mtest <- rbindlist(listP) %>% arrange(list_index)
saveRDS(mtest, "NPEmtest.rds")
# if p > 0.05ish, then no difference between test of null M vs M+1, therefore simpler model is favored.

