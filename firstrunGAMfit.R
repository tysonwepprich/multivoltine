
# file needs to be edited depending on whether Windows or Linux for package loading
source('bootstrapMfunctions.R')



firstrun <- readRDS("ssskipdata.rds")
specBrood <- readRDS("firstrunbrood.rds")

i <- 6
pred <- firstrun[[i]]$GAMpred
species <- specBrood$species[i]
minBrood <- specBrood$minBrood[i]
maxBrood <- specBrood$maxBrood[i]
# species <- "Silver-spotted Skipper"
# minBrood <- 2
# maxBrood <- 3

yr <- as.character(1998:2014)
obsM <- c(minBrood:maxBrood) 
nSim <- as.character(1:5)

params <- expand.grid(species, yr, obsM, nSim,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "yr", "obsM", "nSim")

# data_file Rdata
dataIN8 <- c("pred", "params")
save(list = dataIN8, file = "dataIN8.RData")

# simple param file for slurm.apply
paramIN8 <- data.frame(nRun = seq(1:nrow(params)))

# # single core
# system.time({
#   test <- lapply(paramIN$nRun, StopoverGAM)
# })

# 
# # multicore
# system.time({
#   cl <- makeCluster(8)
#   clusterEvalQ(cl, {
#     library(devtools)
#     library(msm)
#     library(reshape)
#     library(dplyr)
#     # library(StopoverCode) #on linux
#     devtools::load_all("StopoverCode", recompile = TRUE) # on windows
#     load("dataIN.RData")
#   })
#   test <- parLapply(cl, paramIN$nRun, StopoverGAM)
#   stopCluster(cl)
# })


job8 <- slurm_apply(f = StopoverGAM, params = paramIN8, 
                         cpus_per_node = 8, nodes = 5, 
                         data_file = "dataIN8.RData", 
                         # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                         output = "raw")




# extract data from SlurmCov results
slurm_codes <- c("slr6464")
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



# get original params files for the species to guide extracting data from slurm output
limit <- length(slurm_out) + 1
i <- 1
reslist <- list()
while(i < limit){
  j <- length(reslist) + 1
  reslist[[j]] <- NA # initialize?
  if(i < (limit - 1)){
    if(length(slurm_out[[i]]) == 4 & length(slurm_out[[i+1]]) == 18){ 
      reslist[[j]]$pars <- slurm_out[[i]]
      reslist[[j]]$stopoverfit <- slurm_out[[i+1]]
      reslist[[j]]$testpars <- params[j, ]
      i <- i+2
    }else{
      if(length(slurm_out[[i]]) == 1){
        reslist[[j]]$pars <- NA
        reslist[[j]]$stopoverfit <- NA
        reslist[[j]]$testpars <- params[j, ]
        i <- i+1
      }else{
        if(length(slurm_out[[i]]) == 2){
          reslist[[j]]$pars <- slurm_out[[i]]$pars
          reslist[[j]]$stopoverfit <- slurm_out[[i]]$stopoverfit
          reslist[[j]]$testpars <- params[j, ]
          i <- i+1
        }
      }
    }
  }else{
    if(length(slurm_out[[i]]) == 1){
      reslist[[j]]$pars <- NA
      reslist[[j]]$stopoverfit <- NA
      reslist[[j]]$testpars <- params[j, ]
      i <- i+1
    }else{
      if(length(slurm_out[[i]]) == 2){
        reslist[[j]]$pars <- slurm_out[[i]]$pars
        reslist[[j]]$stopoverfit <- slurm_out[[i]]$stopoverfit
        reslist[[j]]$testpars <- params[j, ]
        i <- i+1
      }
    }
  }
}

# now each row of params corresponds to a 2-part list in reslist
outlist <- list()
for (i in 1:length(reslist)){
  templist <- reslist[[i]]
  
  if(is.na(templist$stopoverfit[1])){
    
    pars <- templist$testpars
    outlist[[i]] <- data.frame(species = pars$species,
                               year = pars$yr,
                               obsM = pars$obsM,
                               nSim = pars$nSim,
                               loglik = NA,
                               npar = NA,
                               phi.est = NA,
                               d1.est = NA,
                               b1.est = NA,
                               SiteID = NA,
                               N.est = NA,
                               brood = NA,
                               mu = NA,
                               sigma = NA,
                               w = NA)
    
    
  }else{
  pars <- templist$pars
  loglik <- templist$stopoverfit$ll.val
  npar <- templist$stopoverfit$npar

  sitedf <- data.frame()
  for (j in 1:pars$obsM){
  sitedf <- rbind(sitedf, data.frame(SiteID = names(templist$stopoverfit$N.est),
                       N.est = templist$stopoverfit$N.est,
                       brood = j,
                       mu = templist$stopoverfit$mu.est[, j],
                       sigma = templist$stopoverfit$sigma.est[, j],
                       w = templist$stopoverfit$w.est[, j]
                       ))
  }
  yrdf <- data.frame(species = pars$species,
                     year = pars$yr,
                     obsM = pars$obsM,
                     nSim = pars$nSim,
                     loglik = loglik,
                     npar = npar,
                     phi.est = templist$stopoverfit$phi.est[1],
                     d1.est = templist$stopoverfit$d1.est,
                     b1.est = templist$stopoverfit$b1.est)
  
  outdf <- cbind(yrdf, sitedf)
  outlist[[i]] <- outdf
  }
}

specresults <- data.table::rbindlist(outlist)

a <- specresults %>% 
  group_by(year, obsM) %>%
  summarise(llmean = mean(loglik, na.rm = TRUE),
            ll.na = length(which(is.na(mu))),
            weirdmu = length(which(mu > 50)),
            weirdsigma = length(which(sigma > 10))) %>%
  arrange(year, -llmean)

b <- specresults %>%
  filter(year == 2014, nSim == 1)
index <- b[which(b$loglik == max(b$loglik, na.rm = TRUE))[1], c("year", "obsM", "nSim"), with = FALSE]
plotindex <- which(params$yr == index$year & params$obsM == index$obsM & params$nSim == index$nSim)  
  temp <- reslist[[plotindex]]
  matplot(t(temp$stopoverfit$betta.est), type = "l")
  

# plotindex <- which(b$year == 2014 & b$ll.rank == 1)
plotindex <- which(params$yr == 2014 & params$nSim == 5)
for (i in 1:length(plotindex)){
  temp <- reslist[[plotindex[i]]]
  matplot(t(temp$stopoverfit$betta.est), type = "l")
}


c <- specresults %>%
  group_by(year, obsM, brood, SiteID) %>%
  summarise(meanmu = mean(mu, na.rm = TRUE),
            meansigma = mean(sigma, na.rm = TRUE),
            meanw = mean(w, na.rm = TRUE),
            meanN = mean(N.est, na.rm = TRUE))

# function to fit stopover model to GAM predictions
# Year is character, obsM best guess from GAM plots
StopoverGAM <- function(nRun){
  
  pars <- params[nRun, ]
  sp <- pars$species
  nsim <- pars$nSim
  obsM <- pars$obsM
  yr <- pars$yr
  outlist <- list()
  # for (i in 1:nsim){
    counts <- pred %>%
      filter(Year == yr) %>%
      data.frame()
    
    # if modeled by gdd
    counts <- counts[which(counts$cumdegday %% 50 == 0), ]
    counts$Week <- counts$cumdegday / 50 + 1
    
    counts_uniq <- counts %>% 
      group_by(SiteID, Week) %>%
      arrange(Week) %>%
      summarise(Total = rpois(1, GAM.pred[1])) # here's the random part in each sim
    
    count_matrix <- as.matrix(reshape::cast(counts_uniq, SiteID ~ Week, value = "Total"))
    
    site_covs <- counts %>%
      select(SiteID, lat) %>%
      distinct()
    
    counts <- count_matrix
    
    # dummy p covariate to make model work, but really just want common Np index
    cov.p <- array(rnorm(dim(counts)[1]*dim(counts)[2]), dim = c(dim(counts), 1))
    cov.p[,,1][is.na(counts)] <- NA
    # select sites with enough individuals counted
    siteRows <- which(rowSums(counts, na.rm = TRUE) >= 5)
    counts <- counts[siteRows, ]
    counts[is.na(counts)] <- -1
    
    M <<- obsM
    S <<- dim(counts)[1]
    K <<- dim(counts)[2]
    TIME <<- dim(counts)[2]
    
    cov.p <- cov.p[siteRows, , , drop = FALSE]
    # qp is the number of covariates for p
    qp <- dim(cov.p)[3]
    for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
    cov.p[is.na(cov.p)] <- -1
    
    
    # Time covariate for phi
    cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))),S,K-1,byrow=TRUE) 
    
    # Covariate (latitude)
    # lat <- rowMeans(cov_array[,,5], na.rm = TRUE)
    lat <- site_covs[, "lat"]
    lat <- lat[siteRows]
    cov.w <- cov.mu <- scale(lat)[,1]
    
    # don't know why double assignment makes things work in parallel
    counts <<- counts
    cov.w <<- cov.w
    cov.mu <<- cov.mu
    cov.p <<- cov.p
    cov.phi <<- cov.phi
    qp <<- dim(cov.p)[3]
    
    
    # p.m <- "cov"
    p.m <<- "common"
    w.m <<- "cov"
    mu.m <<- "cov"
    sigma.m <<- "het"
    # phi.m <- "logit.a"
    phi.m <<- "const"
    if (M == 1){ # this needed because code gave errors otherwise
      sigma.m <<- "hom"
      w.m <<- "common"
    } 
    
    
    
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
      tempchoose <- min(c(1:Tries)[which(temp.ll==max(temp.ll, na.rm=TRUE))])
      temp <-temp.fit[[tempchoose]]
    }else{
      temp <- list()
      temp$ll.val <- NA
    }
    
    temp <-temp.fit[[tempchoose]] 
    outlist$pars <- pars
    outlist$stopoverfit <- temp
    
# }
  return(outlist)
}
