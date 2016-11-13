#once new GAMS done, rerun this
#check mclust of sites, changed modelnames for region9 to keep it at 9 with new sites.


rm(list=ls())
library(mclust)
library(plyr)
library(dplyr)
library(mgcv)
library(lubridate)
library(tidyr)
library(stringr)
library(ggplot2)
library(geosphere)
library(data.table)
library(MASS)
library(rslurm)

sites <- read.csv("data/OHsites_reconciled.csv")
names(sites)[1] <- "site"
siteGDD <- readRDS("data/growingDD_Daymet.RDS")
siteGDD <- siteGDD %>%
  group_by(site) %>% 
  filter(yday == 365) %>%
  summarise(meanGDD = mean(cumdegday))
sites <- merge(sites, siteGDD, by = "site")
sitemod <- densityMclust(scale(sites[,c(3:5)]), G = 1:15)
sites$region9 <- as.character(sitemod$classification)
sitemod <- densityMclust(scale(sites[,c(3:5)]), G = 4)
sites$region4 <- as.character(sitemod$classification)

sites$region4 <- plyr::mapvalues(sites$region4, from = c("1", "2", "3", "4"), 
                                 to = c("NE", "NW", "CN", "SW"))
sites$region9 <- plyr::mapvalues(sites$region9, from = c(as.character(1:9)),
                                 to = c("Cuyahoga", "Painesville", "Huron",
                                        "Findlay", "Columbus", "Cincinnati",
                                        "Toledo", "Dayton", "Lancaster"))
sites$SiteID <- formatC(sites$site, width = 3, format = "d", flag = "0")

# plot(sitemod, what = "density")

# visualize clusters
# sites$class <- as.character(sitemod$classification)
# a <- ggplot(data = sites, aes(x = lon, y = lat, group = region9, color = region9)) + geom_point()
# a

site_geo <- read.csv("data/OHsites_reconciled.csv", header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))


# visualize GAM predictions from species models
# preddat <- pred %>% filter(Year == 2003)
# a <- ggplot(data = preddat, aes(x = cumdegday, y = GAM.pred.reg9yr, group = SiteID, color = SiteID))
# a + geom_line() + facet_wrap(~ Reg9Year)

gdd <- readRDS("data/growingDD_Daymet.RDS")
gdd <- gdd[gdd$year >= 1995, ]
gdd <- as.data.frame(gdd)
gdd$date <- strptime(paste(gdd$year, gdd$yday, sep = " "), format = "%Y %j")
gdd$SiteDate <- as.Date(gdd$date)
gdd$date <- NULL

gdd <- gdd %>%
  mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))
gdd <- merge(gdd, site_geo, by = "SiteID")
gdd$SiteYear <- paste(gdd$SiteID, gdd$year, sep = "_")
# gdd_tomerge <- gdd %>% select(SiteYear, yday, cumdegday) %>%
#   dplyr::rename(Ordinal = yday)
# #what about expected gdd by photoperiod/site?
# gdd_left <- gdd %>% 
#   group_by(SiteID, year) %>% 
#   mutate(siteyrtotalgdd = max(cumdegday)) %>% 
#   ungroup() %>% 
#   mutate(actualgddleft = siteyrtotalgdd - cumdegday) %>% 
#   group_by(SiteID, yday) %>% 
#   mutate(expgddleft = mean(actualgddleft)) %>% 
#   ungroup() %>% 
#   mutate(gddmismatch = actualgddleft - expgddleft) %>% 
#   select(SiteYear, yday, cumdegday, actualgddleft, expgddleft, gddmismatch) %>% 
#   dplyr::rename(Ordinal = yday)
# 

# rm(gdd)

# i <- 67
fs <- list.files("gamGDDordinal")
#only try a few species first for mixture model brood separation
mvspec <- c(10, 23, 26, 40, 28, 44, 46, 52, 55, 58, 62, 74, 79, 71, 77, 65)
mvbroodmin <- c(2, 1, 2, 3, 1, 2, 1, 1, 1, 2, 2, 2, 3, 2, 2, 2)
mvbroodmax <- c(3, 2, 3, 4, 2, 3, 2, 2, 2, 3, 3, 3, 4, 3, 3, 3)
specorder <- rank(mvspec)
specdf <- data.frame(mvspec, mvbroodmin, mvbroodmax, specorder)
specdf <- specdf %>% arrange(specorder)

mvbroodmin <- specdf$mvbroodmin
mvbroodmax <- specdf$mvbroodmax
fs <- list.files("../../../nfs/insectmodels-data")



# Estimates brood clusters for each SiteYear, based on simulations from GAM predictions
# changed this to include days outside monitoring to try 
# and ensure broods near boundary counted by mclust
for (i in 1:length(fs)){
# for(i in 1:4){
  # BroodMixMod <- function(i){
  # fs <- list.files("../../../nfs/insectmodels-data/")
  gamlist <- readRDS(paste("../../../nfs/insectmodels-data/", fs[i], sep = ""))
  gammod <- gamlist$modb
  datGAM <- gamlist$counts
  species <- gamlist[[1]]
  rm(gamlist)
  
  mvmin <- mvbroodmin[i]
  mvmax <- mvbroodmax[i]
  
  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday, region9) %>%
    filter(SiteYear %in% unique(datGAM$SiteYear)) %>% 
    # filter(yday >= 75 & yday <= 320) %>%
    filter(yday %in% (seq(62, 328, 7) + sample.int(n=6, size=39, replace=TRUE))) %>%
    dplyr::rename(Ordinal = yday)
  pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year")]))
  
  # pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
  pred <- as.data.frame(pred)
  
  Xp <- predict(gammod, pred, type="lpmatrix") ## map coefs to fitted curves
  beta <- coef(gammod)
  Vb   <- vcov(gammod) ## posterior mean and cov of coefs
  n <- 500 # choose number of simulations
  # set.seed(10)
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  ilink <- family(gammod)$linkinv
  randpred <- array(data = NA, dim = c(n, dim(pred)[1]))
  for (k in seq_len(n)) { 
    modpred   <- ilink(Xp %*% mrand[k, ])
    for (j in seq_len(length(modpred))){
      randpred[k, j] <- rnbinom(1, mu = modpred[j], size = gammod$family$getTheta(TRUE))
    }
  }
  
  paramIN <- data.frame(nsim = c(1:n))
  
  #regional
  SimFunc <- function(nsim){
    ns <- nsim
    outlist <- list()
    for (j in 1:length(unique(pred$region9))){
      index <- unique(pred$region9)[j]
      rowind <- which(pred$region9 == index)
      clustdat <- pred %>%
        filter(region9 == index) %>%
        mutate(GAMpois = randpred[ns, c(rowind)])
      
      if(sum(clustdat$GAMpois) < 10){
        next
      }
      
      dd <- rep(clustdat$cumdegday, clustdat$GAMpois)
      daycount <- rep(clustdat$Ordinal, clustdat$GAMpois)
      SiteYear <- rep(clustdat$SiteYear, clustdat$GAMpois)
      Reg9Year <- rep(clustdat$Reg9Year, clustdat$GAMpois)
      dat <- data.frame(SiteYear, Reg9Year, dd, daycount)

      mod.dd <- try(Mclust(dat$dd, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
      # mod.day <- try(Mclust(dat$daycount, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
      if(class(mod.dd) == "try-error") next
      dat$class <- as.character(mod.dd$classification)
      
      df1 <- data.frame()
      for (b in 1:length(mod.dd$parameters$mean)){
        brooddat <- dat %>% filter(class==b)
        for (c in 1:length(unique(brooddat$Reg9Year))){
          broodtemp <- brooddat %>% filter(Reg9Year == unique(brooddat$Reg9Year)[c])
          if(length(unique(broodtemp$dd))==1){  #mclust can't fit if only one unique dd value
            df1 <- rbind(df1, data.frame(model = "degday",
                                         nsim = ns,
                                         Reg9Year = unique(brooddat$Reg9Year)[c],
                                         brood = b,
                                         num = length(broodtemp$dd),
                                         mu = broodtemp$dd[1],
                                         sigma = NA))
          }else{
            broodmod <- try(Mclust(broodtemp$dd, G=1), silent = TRUE)
            if(class(broodmod)=="try-error") next
            df1 <- rbind(df1, data.frame(model = "degday",
                                         nsim = ns,
                                         Reg9Year = unique(brooddat$Reg9Year)[c],
                                         brood = b,
                                         num = nrow(broodtemp),
                                         mu = broodmod$parameters$mean,
                                         sigma = sqrt(broodmod$parameters$variance$sigmasq)))
          }
        }
      }
      # df2 <- data.frame()
      # for (c in 1:length(mod.day$parameters$mean)){
      #   df2 <- rbind(df2, data.frame(model = "ordinal",
      #                                nsim = ns,
      #                                SiteYear = index,
      #                                brood = c,
      #                                num = length(which(mod.day$classification == c)),
      #                                weight = mod.day$parameters$pro[c],
      #                                mu = mod.day$parameters$mean[c],
      #                                sigma = sqrt(mod.day$parameters$variance$sigmasq[1]))
      #   )
      # }
      # df <- rbind(df1,df2)
      outlist[[length(outlist) + 1]] <- df1
    }
    outdf <- data.table::rbindlist(outlist)
    
    pops <- clustdat %>% 
      group_by(SiteYear) %>% 
      summarise(SiteN=sum(GAMpois)) %>% 
      mutate(nsim=ns)
    
    outfunc <- list("brood" = outdf, "popsize" = pops)
    
    return(outfunc)
  }
  # original
  # SimFunc <- function(nsim){
  #   ns <- nsim
  #   outlist <- list()
  #   for (j in 1:length(unique(pred$SiteYear))){
  #     index <- unique(pred$SiteYear)[j]
  #     clustdat <- pred %>%
  #       filter(SiteYear == index) %>%
  #       # select(Ordinal, GAM.pred, lat, SiteID) %>%
  #       dplyr::select(Ordinal, cumdegday, GAM.pred.reg9yr) %>%
  #       rowwise() %>%
  #       mutate(GAMpois = rpois(1, GAM.pred.reg9yr[1])) %>%
  #       ungroup()
  #     
  #     if(sum(clustdat$GAMpois) < 10){
  #       next
  #     }
  #     
  #     dd <- rep(clustdat$cumdegday, clustdat$GAMpois)
  #     daycount <- rep(clustdat$Ordinal, clustdat$GAMpois)
  #     dat <- data.frame(dd, daycount)
  #     # 
  #     # t <-  try(modbs <- mclustBootstrapLRT(dat$dd, model = "E",
  #     #                              nboot = 500,
  #     #                              verbose = FALSE))
  #     mod.dd <- try(Mclust(dat$dd, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
  #     mod.day <- try(Mclust(dat$daycount, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
  #     if(class(mod.dd) == "try-error") next
  #     # dat$class <- as.character(mod$classification)
  #     
  #     df1 <- data.frame()
  #     for (b in 1:length(mod.dd$parameters$mean)){
  #       df1 <- rbind(df1, data.frame(model = "degday",
  #                                    nsim = ns,
  #                                    SiteYear = index,
  #                                    brood = b,
  #                                    num = length(which(mod.dd$classification == b)),
  #                                    weight = mod.dd$parameters$pro[b],
  #                                    mu = mod.dd$parameters$mean[b],
  #                                    sigma = sqrt(mod.dd$parameters$variance$sigmasq[1]))
  #       )
  #     }
  #     df2 <- data.frame()
  #     for (c in 1:length(mod.day$parameters$mean)){
  #       df2 <- rbind(df2, data.frame(model = "ordinal",
  #                                    nsim = ns,
  #                                    SiteYear = index,
  #                                    brood = c,
  #                                    num = length(which(mod.day$classification == c)),
  #                                    weight = mod.day$parameters$pro[c],
  #                                    mu = mod.day$parameters$mean[c],
  #                                    sigma = sqrt(mod.day$parameters$variance$sigmasq[1]))
  #       )
  #     }
  #     df <- rbind(df1,df2)
  #     outlist[[length(outlist) + 1]] <- df
  #   }
  #   outdf <- data.table::rbindlist(outlist)
  #   return(outdf)
  # }
  
  ty <- slurm_apply(f = SimFunc, params = paramIN, nodes = 4,
                    cpus_per_node = 8, 
                    add_objects = c("pred", "mvmin", "mvmax", "randpred"),
                    slurm_options = list(partition = "sesync"))
}

fs <- list.files("../../../nfs/insectmodels-data/")
for (i in 1:length(fs)){
  gamlist <- readRDS(paste("../../../nfs/insectmodels-data/", fs[i], sep = ""))
  slurm_call(f = BroodMixMod, 
             add_objects = c("gdd", "mvbroodmax", "mvbroodmin", "gamlist"))
}

library(parallel)
# multicore
system.time({
  cl <- makeCluster(4)
  clusterEvalQ(cl, {
    library(mclust)
    library(data.table)
    library(dplyr)
    library(mgcv)
  })
  clusterExport(cl=cl,
                varlist=c("mvspec",
                          "mvbroodmax",
                          "mvbroodmin",
                          "gdd"))
  
  test <- parLapply(cl, 1:16, BroodMixMod)
  stopCluster(cl)
})



# pulling  results together
resdf <- data.frame()
resfile <- list.files()[grep("_rslurm_slr", list.files())]
for (i in 1:length(resfile)){
  tempdir <- resfile[i]
  # jobnum <- as.numeric(gsub("_rslurm_slr", "", tempdir))
  timecreate <- file.info(paste(tempdir, "submit.sh", sep = "/"))$ctime
  resdf <- rbind(resdf, data.frame(tempdir, timecreate))
}

resdf <- resdf %>% arrange(timecreate)
results <- cbind(specdf, resdf[59:74,])
specname <- gsub("gam.gdd.plus.day", "", fs)
CommonName <- gsub(".rds", "", specname)
results$CommonName <- CommonName

allresults <- list()
for (i in 1:16){
  for(j in 0:15){
    resname <- paste("results", j, sep = "_")
    resname <- paste(resname, ".RData", sep = "")
    try(temp <- readRDS(paste(results$tempdir[i], resname, sep = "/")))
    if (class(temp) == "character") next
    allbroods <- list()
    allpops <- list()
    for(k in 1:length(temp)){
      allbroods[[k]] <- temp[[k]]$brood
      allpops[[k]] <- temp[[k]]$pop
    }
    outdf1 <- rbindlist(allbroods)
    outdf2 <- rbindlist(allpops)
    outdf1$species <- results$CommonName[i]
    outdf2$species <- results$CommonName[i]
    allresults[[length(allresults) + 1]] <- list("brood"=outdf1, "pop"=outdf2)
  }
}
# alldf <- rbindlist(allresults)
# saveRDS(alldf, "mclust_results_reg9yr0.rds")
saveRDS(allresults, "mclust_results_region9.rds")

# can't read this rds on my laptop
test <- readRDS("mclust_results_region9.rds") #weird duplicates in the later nsims???

brood <- list()
pop <- list()
for (i in 1:length(test)){
  brood[[i]] <- test[[i]]$brood
  pop[[i]] <- test[[i]]$pop
}
brooddf <- data.table::rbindlist(brood)
popdf <- data.table::rbindlist(pop)
brooddf <- dplyr::distinct(brooddf)
popdf <- dplyr::distinct(popdf)
saveRDS(brooddf, file="mclust_broods.rds")
saveRDS(popdf, file="mclust_pops.rds")
#### one issue with regional clustering, no site population estimates


# Estimates brood clusters for each SiteYear, based on simulations from GAM predictions
# changed this to include days outside monitoring to try 
# and ensure broods near boundary counted by mclust
specieslist <- list()
for (i in 1:length(fs)){
  # BroodMixMod <- function(i){
  # fs <- list.files("../../../nfs/insectmodels-data/")
  gamlist <- readRDS(paste("../../../nfs/insectmodels-data/", fs[i], sep = ""))
  gammod <- gamlist$modb
  datGAM <- gamlist$counts
  species <- gamlist[[1]]
  rm(gamlist)
  
  mvmin <- mvbroodmin[i]
  mvmax <- mvbroodmax[i]
  
  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday) %>%
    filter(SiteYear %in% unique(datGAM$SiteYear)) %>% 
    filter(yday >= 75 & yday <= 320) %>%
    # filter(yday %in% seq(55, 335, 7)) %>% 
    dplyr::rename(Ordinal = yday)
  pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year")]))
  
  # pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
  pred <- as.data.frame(pred)
  
  Xp <- predict(gammod, pred, type="lpmatrix") ## map coefs to fitted curves
  beta <- coef(gammod)
  Vb   <- vcov(gammod) ## posterior mean and cov of coefs
  n <- 100
  set.seed(10)
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  ilink <- family(gammod)$linkinv
  randpred <- array(data = NA, dim = c(n, dim(pred)[1]))
  for (k in seq_len(n)) { 
    modpred   <- ilink(Xp %*% mrand[k, ])
    for (j in seq_len(length(modpred))){
      randpred[k, j] <- rnbinom(1, mu = modpred[j], size = gammod$family$getTheta(TRUE))
    }
  }
  
  #regional
  outlist <- list()
  
  for (ns in 1:n){    
    for (j in 1:length(unique(pred$Reg9Year))){
      index <- unique(pred$Reg9Year)[j]
      rowind <- which(pred$Reg9Year == index)
      clustdat <- pred %>%
        filter(Reg9Year == index) %>%
        mutate(GAMpois = randpred[ns, c(rowind)]) %>% 
        group_by(SiteYear) %>% 
        summarise(totalN = sum(GAMpois))
      clustdat$Reg9Year <- index
      clustdat$nsim <- ns
      clustdat$species <- species
      outlist[[length(outlist) + 1]] <- clustdat
    }
  }
  outdf <- data.table::rbindlist(outlist)
  specieslist[[length(specieslist) + 1]] <- outdf
}
finaldf <- data.table::rbindlist(specieslist)

pops <- finaldf %>% 
  group_by(species, SiteYear) %>% 
  summarise(meanN = mean(totalN),
            sdN = sd(totalN),
            medN = median(totalN))
saveRDS(pops, "provisionalN_SiteYear.rds")
