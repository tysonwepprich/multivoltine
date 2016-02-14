# source this for bootstrap M functions
list.of.packages <- c("devtools", "msm", "parallel", "plyr", "dplyr", "tidyr", 
                      "readr", "reshape", "reshape2", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(rslurm)
library(devtools)
library(msm)
library(parallel)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(reshape)
library(reshape2)
library(data.table)


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


# Data prep function, input CommonName
# Output list of counts, survey covariates, and site covariates for each year 1998-2012
SpeciesData <- function(species){
  # Array of species counts (1 species x sites x years)
  # Covariate array to match
  
  ##########
  #DATA PREP
  ##########
  
  data <- fread("data/data.trim.csv", header = TRUE)
  # data <- data[, list(SeqID, SiteID.x, SiteDate, Week, Total, CheckListKey, CommonName)]
  setnames(data,"SiteID.x","SiteID")
  data[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]
  data[, SiteDate := parse_date(SiteDate, format = "%Y-%m-%d")]
  
  #Site covariates for mu, weights
  #GDD from Dan (or maybe Leslie/Rick Reeves?)
  gdd <- fread("data/GddResultsAllSites_1996_2012.csv")
  names(gdd)[1] <- "SiteID"
  gdd$SiteID <- formatC(gdd$SiteID, width = 3, format = "d", flag = "0")
  gdd$Date <- parse_date(gdd$full_date, format = "%m/%d/%Y")
  gdd$Year <- year(gdd$Date)
  gdd_summary <- gdd %>%
    filter(todayGDD >= 0) %>%
    group_by(SiteID, Year) %>%
    summarise(YearGDD = max(GDD),
              SpringGDD = max(GDD[ordinalEndDayOfYear < 100], na.rm = TRUE), #negative infinity popping up????
              SummerGDD = max(GDD[ordinalEndDayOfYear < 200], na.rm = TRUE),
              FallGDD = max(GDD[ordinalEndDayOfYear < 300], na.rm = TRUE))
  
  gdd_summary[SpringGDD == -Inf]$SpringGDD <- NA
  
  
  sites <- read_csv("data/OHsites_reconciled.csv")
  names(sites)[1] <- "SiteID"
  sites$SiteID <- formatC(sites$SiteID, width = 3, format = "d", flag = "0")
  
  # resolve sites not matching between 3 datasets
  masterSites <- merge(unique(sites[, "SiteID"]), unique(gdd_summary[, "SiteID", with = FALSE]))
  masterSites <- merge(masterSites, unique(data[, "SiteID", with = FALSE]))
  
  data <- data[which(data$SiteID %in% masterSites$SiteID), ]
  sites <- sites[which(sites$SiteID %in% masterSites$SiteID), ]
  gdd_summary <- gdd_summary[which(gdd_summary$SiteID %in% masterSites$SiteID), ]
  
  data[, Year := year(SiteDate)]
  data[, `:=` (WeekPerYear = length(unique(Week)),
               SurvPerYear = length(unique(SeqID))), 
       by = list(SiteID, Year)]
  dat <- data[WeekPerYear >= 15]
  
  surveys <- distinct(dat[, c("SiteID", "SiteDate", "Week", "SeqID", "Year"), with = FALSE])
  dat <- dat[CommonName == species][Year >= 1998]
  
  
  years <- sort(unique(dat$Year))
  dat_list <- as.list(years)
  for (i in 1:length(dat_list)){
    
    yr <- years[i]
    spdat <- dat[Year == yr]
    spdat <- unique(spdat)
    #   setkey(spdat, SeqID)
    
    survs <- surveys[year(SiteDate) == yr]
    #Add zeros to surveys when species not counted during a survey
    
    test <- merge(survs, spdat, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"), all.x = TRUE)
    counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year"), with = FALSE]
    counts$Total <- mapvalues(counts[,Total], from = NA, to = 0)
    
    #overachieving volunteers going out more than once a week!
    #choose first one, not averaging, easier to match with covariates
    counts_uniq <- counts %>% 
      group_by(SiteID, Week) %>%
      arrange(SiteDate) %>%
      summarise(Total = Total[1])
    
    count_matrix <- as.matrix(cast(counts_uniq, SiteID ~ Week, value = "Total"))
    count_matrix <- round(count_matrix)
    count_matrix[is.na(count_matrix)] <- NA
    
    dat_list[[i]]$counts <- count_matrix
    
    
    #covariates
    
    #some already calculated in OHdetprob.RMD
    #what to do about NA's in covariates?
    oldcovs <- fread("data/survey.covariates.csv")
    covs <- merge(survs, oldcovs, by = "SeqID", all.x = TRUE)
    
    #Celsius-Fahrenheit issues
    covs[mean.temp < 45]$mean.temp <- covs[mean.temp < 45]$mean.temp * 1.8 + 32
    
    #overachieving volunteers going out more than once a week!
    #choose first one, not averaging, to match with counts
    covs <- covs %>% 
      group_by(SiteID, Week) %>%
      arrange(SiteDate) %>%
      mutate(Duplicate = 1:length(SeqID)) %>%
      filter(Duplicate == 1)
    
    surv <- survs %>%
      group_by(SiteID, Week) %>%
      arrange(SiteDate) %>%
      summarise(SiteDate = SiteDate[1])
    
    covs <- merge(surv, covs, by = c("SiteID", "Week", "SiteDate"), all.x = TRUE)
    
    # some NA's, not more than 30 for covariates
    # just assign them as mean (even though not perfect for seasonal variables)
    covs <- data.frame(covs)
    for(j in 6:ncol(covs)){
      covs[is.na(covs[,j]), j] <- mean(covs[,j], na.rm = TRUE)
    }
    
    covs$Ztemp <- scale(poly(covs$mean.temp, 2)[, 1])[1:nrow(covs)]
    covs$Ztemp2 <- scale(poly(covs$mean.temp, 2)[, 2])[1:nrow(covs)]
    covs$Zwind <- scale(covs$mean.wind)[1:nrow(covs)]
    covs$Zcloud <- scale(covs$mean.cloud)[1:nrow(covs)]
    covs$Zduration <- scale(covs$duration)[1:nrow(covs)]
    covs$Zhour <- scale(covs$start.hour)[1:nrow(covs)]
    covs$Zspecies <- scale(covs$num.species)[1:nrow(covs)]
    covs$Zallabund <- scale(log(covs$abund + 1))[1:nrow(covs)]
    
    # new idea for phi, include ordinal day for time/age standard instead of week
    covs$Zjulian <- scale(yday(covs$SiteDate))
    
    #cast covs as matrix, so NA's inserted for missing surveys
    cov_array <- array(NA, dim=c(length(unique(surv$SiteID)), length(unique(surv$Week)), 9))
    
    cov_molten <- melt(covs, id = c("SiteID", "Week", "SiteDate"))
    cov_array[,,1] <- as.matrix(cast(cov_molten[cov_molten$variable == "Ztemp", ], SiteID ~ Week, value = "value"))
    cov_array[,,2] <- as.matrix(cast(cov_molten[cov_molten$variable == "Ztemp2", ], SiteID ~ Week, value = "value"))
    cov_array[,,3] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zwind", ], SiteID ~ Week, value = "value"))
    cov_array[,,4] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zcloud", ], SiteID ~ Week, value = "value"))
    cov_array[,,5] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zduration", ], SiteID ~ Week, value = "value"))
    cov_array[,,6] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zhour", ], SiteID ~ Week, value = "value"))
    cov_array[,,7] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zspecies", ], SiteID ~ Week, value = "value"))
    cov_array[,,8] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zjulian", ], SiteID ~ Week, value = "value"))
    cov_array[,,9] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zallabund", ], SiteID ~ Week, value = "value"))
    
    dat_list[[i]]$surv_covs <- cov_array
    

    
    #spring and yearly GDD have lowest correlation, but still .62.
    gdd_covs <- gdd_summary[Year == yr]
    
    cov_sites <- merge(sites, gdd_covs, by = "SiteID", all.x = TRUE)
    
    # gdd missing from Catawba Island 103
    # use gdd from closest other site
    # Turn this into a function!
    rowNA <- which(is.na(cov_sites$SpringGDD))
    if (length(rowNA) > 0){
      d <- dist(cbind(cov_sites$lat, cov_sites$lon), upper = TRUE)
      dists <- as.matrix(d)[rowNA,]
      dists[which(dists == 0)] <- NA
      mindist <- apply(dists, 1, min, na.rm = TRUE)
      for (md in 1:length(rowNA)){
        rowReplace <- which(dists[md, ] == mindist[md])
        cov_sites[rowNA[md], c("Year","YearGDD", "SpringGDD", "SummerGDD", "FallGDD")] <-  cov_sites[rowReplace, c("Year", "YearGDD", "SpringGDD", "SummerGDD", "FallGDD")]
      }
    }
    cov_sites <- cov_sites[which(cov_sites$SiteID %in% unique(survs$SiteID)), ]
    
    dat_list[[i]]$site_covs <- cov_sites
    
  }
  return(dat_list)
}

# slurm.apply function
# fits each set of covariates once (multiple tries to find global LL)

SlurmCovs <- function(nRun){
  # parameters
  # is <<- necessary so inner functions can find all these?
  pars <- params[nRun, ]
  species <- pars$species
  list_index <- pars$list_index
  raw_cutoff <- pars$raw_cutoff
  p_cov1 <- pars$p_cov1
  p_cov2 <- pars$p_cov2
  p_cov3 <- pars$p_cov3
  p_cov4 <- pars$p_cov4
  site_covs <- pars$site_covs
  M <<- pars$M
  #   p.m <<- pars$p.m
  #   w.m <<- pars$w.m
  #   mu.m <<- pars$mu.m
  sigma.m <<- pars$sigma.m
  phi.m <<- pars$phi.m
  
  ### data prep ###
  temp <- dat[[list_index]]
  year <- temp[[1]]
  counts <- temp$counts
  
  # select sites with enough individuals counted
  siteRows <- which(rowSums(counts, na.rm = TRUE) >= raw_cutoff)
  surv_present <- which(apply(counts, 1, function(x) length(which(x > 0))) >= 3)
  siteRows <- siteRows[which(siteRows %in% surv_present)]
  counts <- counts[siteRows, ]
  counts[is.na(counts)] <- -1 #reassign NA for C program
  ALLcounts <- counts # ALL distinguishes from bootstrapped counts
  
  S <<- dim(counts)[1] 
  K <<- dim(counts)[2]  
  TIME <<- dim(counts)[2]
  
  # covariates
  cov_array <- temp$surv_covs
  # detection probability
  covs <- c(p_cov1, p_cov2, p_cov3, p_cov4)
  covs <- sort(as.numeric(covs[covs %in% as.character(c(1:9))]))
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
  
  cov_sites <- temp$site_covs
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
  temp$nRun <- nRun
  out[[1]] <- temp
  return(out)
}


# Fit models for null and alternative (M + 1) for each simulated dataset
SlurmGeneration <- function(nRun){
  out <- vector("list", length = 2)
  BSdata <- SampleList[[nRun]]
  if (length(which("simData" %in% names(BSdata))) == 0){ # if no simData made, move along (happens if ll.val = NA or rpois error (N.est too big))
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
