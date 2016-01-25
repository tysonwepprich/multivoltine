.tmplib <- lapply(c('base','methods','datasets','utils','grDevices','graphics','stats','rslurm','devtools','msm','parallel','plyr','dplyr','tidyr','readr','reshape','reshape2','data.table','StopoverCode'), 
           library, character.only = TRUE, quietly = TRUE) 
load('dataIN.RData') 
load('slr7310.RData') 
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
.rslurm_istart <- .rslurm_id * 700 + 1 
.rslurm_iend <- min((.rslurm_id + 1) * 700, 
                    nrow(.rslurm_params)) 
.rslurm_result <- do.call(parallel::mcmapply, c(function(nRun){
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
, .rslurm_params[.rslurm_istart:.rslurm_iend, , drop = FALSE], mc.cores = 8)) 
save(.rslurm_result, file = paste0('slr7310_', .rslurm_id, '.RData'))