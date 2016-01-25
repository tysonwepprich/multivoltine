# abbreviated version of bootstrampM.R
# run one M null hypothesis at a time, to prevent unnecessary computation

source('bootstrapMfunctions.R')


# load data
# count_array <- readRDS('count_array_expanded.rds')
# cov_array <- readRDS('covariates_array_expanded.rds')
# cov_sites <- readRDS('covariates_sites_expanded.rds')
# species_list <- readRDS('species_expanded.rds')
# 
count_array <- readRDS('count_array.2010.rds')
cov_array <- readRDS('covariates_array.2010.rds')
cov_sites <- readRDS('covariates_sites.2010.rds')
species_list <- readRDS('species.2010.rds')

# choose parameter ranges
species <- c(7,8,12:15,17:19,21:23,26,27,29,32,36,38,39) # corresponds to row in species_list
raw_cutoff <- 5
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- "AnnGDD" # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(3, 4) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("count_array", "cov_array", "cov_sites", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

# calculate null hypotheses for M for different species
baseline <- slurm_apply(f = SlurmCovs, params = paramIN, 
                        cpus_per_node = 8, nodes = 5, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")

slurm_codes <- c("slr1539", "slr4788")
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

outDF <- do.call("rbind", outDF)
baselineDF <- outDF
# saveRDS(baselineDF, file = "baseline3v4.rds")


# select M=2 to simulate data (to test vs. M = 3)
# for M 3v4, remove SSSkip, too long to run
index <- which(baselineDF$M == 3 & baselineDF$species != 38) # Northern Broken Dash N.est too high, rpois errors
# index <- which(baselineDF$M == 3)

slurm_out2 <- slurm_out[c(index)]



# what to do about ll.val == NA for M = 1?
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

# saveRDS(SampleList, file = "SampleList.rds")



# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:length(SampleList)))



# calculate null hypotheses for M = 1:5 for different species
alt <- slurm_apply(f = SlurmGeneration, params = paramIN, 
                   cpus_per_node = 8, nodes = 8, 
                   data_file = "dataIN.RData", 
                   # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                   output = "raw")

slurm_codes <- c("slr3887", "slr860")
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


# baselineDF <- readRDS("baseline3v4.rds")
parsIN <- data.frame(expand.grid.alt(c(2,3), unique(BSmods$species)))
names(parsIN) <- c("nullM", "spec")
Mtest <- parsIN %>% rowwise() %>% mutate(pval = BSpval(nullM, spec)) %>% data.frame()
# saveRDS(Mtest, file = "Mtest2v3.rds")
saveRDS(Mtest, file = "Mtest3v4.rds")
