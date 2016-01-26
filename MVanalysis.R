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

# i <- 16 #slr4158
# i <- 15 #slr7296 #Silver Sp Skip has 3M SlurmCov error "Error in rowSums(counts, na.rm = TRUE) : \n  'x' must be an array of at least two dimensions
# i <- 14 #slr7389 #RSP has NA ll.val on 2004/06 for 2/3M, unknown why it's not fitting
# i <- 13 #slr1826
# i <- 12 #slr1965
# i <- 10 #slr2023
# i <- 2 #slr2085
i <- 3 #slr2152

species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i]

# somewhat unwieldly, list with each year as a list of 4 (year, counts, surv_covs, site_covs)
dat <- SpeciesData(species)

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

list_index_min_data <- unique(data_avail$list_index[data_avail$both_met >= 10])

# choose parameter ranges
raw_cutoff <- 5 # c(5, 10)
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- c("AnnGDD", "lat") # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(minBrood:maxBrood) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, list_index_min_data, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "list_index", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("dat", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))


# calculate null hypotheses for M for different species
ETigSwalCovs <- slurm_apply(f = SlurmCovs, params = paramIN, 
                        cpus_per_node = 8, nodes = 2, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")
# 
# cl <- makeCluster(4)
# clusterEvalQ(cl, {
#   library(devtools)
#   library(msm)
#   library(dplyr)
#   library(StopoverCode) #on linux
#   # devtools::load_all("StopoverCode", recompile = TRUE) # on windows
#   load("dataIN.RData")
# })
# test <- parLapply(cl, paramIN$nRun, SlurmCovs)
# stopCluster(cl)

# Next step, go through processSlurmCov.R to choose best site_cov for species over all years
# May need to edit code depending on how lists are output, lengths may differ
# this throws indexing off for some species/years
# Stupidly manual, but at least allows for some checks.




# simulate data from best-fit model parameters for each year
dat <- SpeciesData(species)
nsim <- 100
SampleList <- vector("list", length = length(slurm_out2)*nsim)
# make bootstrap data simulations for each model x 100
for (bs in 1:length(SampleList)){
  mod <- (bs + nsim - 1) %/% nsim        #ADD IN INDEX FOR MOD AND BOOTSTRAP TO TRACK THESE
  nullFit <- slurm_out2[[mod]]
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
  list_index <- nullFit$pars$list_index
  counts <- dat[[list_index]]$counts
  raw_cutoff <- nullFit$pars$raw_cutoff
  
  # select sites with enough individuals counted
  
  siteRows <- which(rowSums(counts, na.rm = TRUE) >= raw_cutoff)
  surv_present <- which(apply(counts, 1, function(x) length(which(x > 0))) >= 3)
  siteRows <- siteRows[which(siteRows %in% surv_present)]
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
  
  # problem with large N.est creating NA's in rpois
  #   if (length(is.na(rpois(S, N.est))) > 0){
  #     SampleList[[bs]]$ll.val <- NA
  #     SampleList[[bs]]$npar <- NA
  #     next
  #   }
  
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

# saveRDS(SampleList, file = "simDataGenMode/SpiceSwalSampleList.rds")
# saveRDS(SampleList, file = "simDataGenMode/PeckSkipSampleList.rds")
saveRDS(SampleList, file = "simDataGenMode/ETigerSwalSampleList.rds")


SampleList <- readRDS("simDataGenMode/ETigerSwalSampleList.rds")
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:length(SampleList)))


cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)



# calculate null hypotheses for same species, different years
peckskip <- slurm_apply(f = SlurmGeneration, params = paramIN, 
                   cpus_per_node = 8, nodes = 2, 
                   data_file = "dataIN.RData", 
                   # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                   output = "raw")


