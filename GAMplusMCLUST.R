#script combining:
#counts filtering
#gam fitting and summary stats
#mclust broods

library(tidyverse)
library(mgcv)
library(rslurm)
library(mclust)
library(data.table)
library(lubridate)
library(MASS)



data <- fread("data/data.trim.csv", header = TRUE)
# data <- data.table(data)
setnames(data,"SiteID.x","SiteID")
data[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]
data$SiteDate <- lubridate::ymd(as.character(data$SiteDate))
data$CommonName <- as.character(data$CommonName)
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week"), with = FALSE])

allspec <- unique(data$CommonName)
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspec[1:122])),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct()

sites <- read.csv("data/OHsites_reconciled_update2016.csv")
names(sites)[1] <- "site"
siteGDD <- readRDS("data/growingDD_Daymet.RDS")
siteGDD <- siteGDD %>%
  group_by(site) %>% 
  filter(yday == 365) %>%
  summarise(meanGDD = mean(cumdegday))
sites <- merge(sites, siteGDD, by = "site")
sitemod <- densityMclust(scale(sites[,c(3:5)]), G = 1:10, modelNames = "VVV")
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

# visualize regional clusters
# sites$class <- as.character(sitemod$classification)
# a <- ggplot(data = sites, aes(x = lon, y = lat, group = region9, color = region9)) + geom_point()
# a

site_geo <- fread("data/OHsites_reconciled_update2016.csv")
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))


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
# gdd_tomerge <- gdd %>% dplyr::select(SiteYear, yday, cumdegday) %>%
#   dplyr::rename(Ordinal = yday)


# 
# SpeciesList <- readRDS("data/SpeciesList.rds")
# species <- SpeciesList %>% arrange(Present) %>% select(CommonName)

speciesdat <- read.csv("data/speciesphenology.csv", header = TRUE)
speciesdat <- speciesdat %>%
  dplyr::select(CommonName, BroodsGAMmin, BroodsGAMmax, UseMV, UseMismatch, SyncedBroods, Notes, Model, OutlierCutoff) %>% 
  filter(UseMismatch == "y")


models <- "extra" #c("orig", "extra", "orig.cov", "extra.cov")
cutoff <- "adapt" #c("adapt","strict", "loose")
params <- expand.grid(speciesdat$CommonName, models, cutoff,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model", "cutoff")
params <- params[c(5),]

# gamlist <- as.list(species$CommonName)
# for (i in 1:nrow(params)){
# i <- 1
# species <- params$species[i]
# model <- params$model[i]
# cutoff <- params$cutoff[i]
FitGAM <- function(species, model, cutoff){
  species <- species
  model <- model
  cutoff <- cutoff
  reduced <- NA
  pars <- data.frame(species, model, cutoff, reduced)
  # pars <- params[nRow,]
  # sp <- pars$species
  # model <- pars$model
  counts <- data %>% filter(CommonName == species) %>% data.table()
  counts[, Ordinal := yday(SiteDate)]
  counts[, Year := year(SiteDate)]
  
  #get unique surveys, including those where species not counted
  survs <- surveys[year(SiteDate) %in% unique(counts$Year)]
  survs <- survs[SiteID %in% unique(counts$SiteID)]
  survs$Year <- year(survs$SiteDate)
  
  #Add zeros to surveys when species not counted during a survey
  test <- merge(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"), all.x = TRUE)
  counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year"), with = FALSE]
  counts$Total <- plyr::mapvalues(counts[,Total], from = NA, to = 0)
  counts <- merge(counts, site_geo, by = "SiteID")
  counts <- merge(counts, covdata, by = "SeqID", all.x = TRUE, all.y = FALSE)
  counts$temperature[which(counts$temperature < 50)] <- NA
  counts$duration[which(counts$duration == 0)] <- NA
  
  #scaling covariates
  # #scaled over whole season/state
  # counts$listlength <- scale(counts$listlength)
  # counts$temperature <- scale(counts$temperature)
  # counts$duration <- scale(counts$duration)
  #scaled by region9/week
  counts <- counts[, `:=` (zlistlength = as.numeric(scale(log(listlength+1)))),
                   by = list(region9, Week)]
  counts <- counts[, `:=` (
    ztemperature = as.numeric(scale(temperature))),
    by = list(region9, Week)]
  # counts <- counts[, `:=` (ztemperature = as.numeric(scale(temperature)))]
  counts <- counts[, `:=` (zduration = as.numeric(scale(duration))),
                   by = list(SiteID)]
  
  
  # trying to add GDD instead of ordinal date
  counts <- merge(counts, gdd, by = c("SiteID", "SiteDate", "lat", "lon", "region9", "region4"), all.x = TRUE, all.y = FALSE)
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  
  
  if(cutoff == "strict"){
    datGAM <- counts[YearTotal >= 3]
    datGAM <- datGAM[SurvPerYear >= 15]
  }
  if(cutoff == "loose"){
    datGAM <- counts[YearTotal >= 1]
    datGAM <- datGAM[SurvPerYear >= 10]
  }
  #to decrease computation of common species, adapt filter removes low data siteyears
  if(cutoff == "adapt"){
    datGAM <- counts[YearTotal >= 1]
    datGAM <- datGAM[SurvPerYear >= 11] 
    #this SurvPerYear cutoff chosen to include 95% of survey effort
    #16 surveys would be including 90% of surveys 
    counts1 <- datGAM %>% 
      group_by(SiteID, Year) %>% 
      mutate(SiteYearSurvs = length(unique(SeqID)))
    counts1 <- counts1 %>% 
      group_by(SiteID) %>% 
      mutate(SpeciesSiteYears = length(unique(Year)))
    counts2 <- counts1 %>% 
      group_by(SiteID, Year, SiteYearSurvs, SpeciesSiteYears) %>% 
      summarise(PresenceSiteYear = length(unique(SeqID[Total>=1])),
                SiteYearTotal = sum(Total))
    
    filtered <- counts2 %>% 
      filter(SpeciesSiteYears >= quantile(counts2$SpeciesSiteYears, 0.1),
             PresenceSiteYear >= quantile(counts2$PresenceSiteYear, 0.1),
             SiteYearTotal >= quantile(counts2$SiteYearTotal, 0.1)) %>% 
      mutate(SiteYear = paste(SiteID, Year, sep = "_"))
    
    datGAM <- datGAM %>% 
      filter(SiteYear %in% unique(filtered$SiteYear))
  }
  
  if(nrow(datGAM) == 0){
    mod <- NA
  }else{
    # datGAM[, SitesObserved := length(unique(SiteID)), by = list(Year)]
    # datGAM <- datGAM[SitesObserved >= 1]
    
    dat <- datGAM
    outlist <- list()
    for (y in sort(unique(dat$Year))){
      temp <- dat %>% filter(Year == y) %>%
        mutate(Ordinal = yday(SiteDate)) %>% 
        group_by(SiteID, Year) %>%
        mutate(YearTotal = sum(Total),
               SurvSeen = length(which(Total > 0))) %>%
        filter(YearTotal > 0,
               SurvSeen > 0) %>%
        dplyr::select(SiteID, Total, zlistlength, lat, lon, region9, Ordinal, 
                      Year, ztemperature, zduration, cumdegday) %>% 
        data.frame()
      
      # pad zeros at beginning and end for GAM fitting
      zeros <- c(30, 40, 50, 60, 70, 320, 330, 340, 350, 360)
      tempdf <- expand.grid(unique(temp$SiteID), zeros)
      names(tempdf) <- c("SiteID", "yday")
      tempdf$SiteID <- as.character(tempdf$SiteID)
      tempdf$year <- y
      tempdf$zlistlength <- 0
      tempdf$ztemperature <- 0
      tempdf$zduration <- 0
      tempdf$Total <- 0
      tempdf <- left_join(tempdf, gdd[,c("SiteID", "year", "yday", "cumdegday")], 
                          by = c("SiteID", "year", "yday"))
      names(tempdf)[2] <- "Ordinal"
      tempdf$year <- NULL
      # tempdf$cumdegday <- c(rep(0, nrow(tempdf)/2), rep(max(temp$cumdegday), nrow(tempdf)/2))
      
      test <- merge(tempdf, distinct(dplyr::select(temp, SiteID, region9, lat, lon, Year)),
                    by = "SiteID", all.x = TRUE, all.y = FALSE)
      
      outlist[[y]] <- rbind(temp, test)
    }
    
    dat <- rbindlist(outlist)
    dat$Year <- as.factor(as.character(dat$Year))
    dat$SiteID <- as.factor(as.character(dat$SiteID))
    dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
    dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
    dat$ztemperature[which(is.na(dat$ztemperature))] <- 0
    dat$zduration[which(is.na(dat$zduration))] <- 0
    dat$Reg9Year <- as.factor(paste(dat$region9, dat$Year, sep = "_"))
    dat <- as.data.frame(dat)
    #silly filters for univoltine species with outliers
    
    if(species == "Baltimore"){
      dat <- dat[-which(dat$cumdegday > 1000 & dat$Total >= 1), ]
    }
    if(species == "Leonard's Skipper"){
      dat <- dat[-which(dat$cumdegday < 1000 & dat$Total >= 1), ]
    }
    
    starttime <- Sys.time()
    
    temp <- dat
    if(sum(temp$Total) < 20|length(unique(temp$SiteID)) < 2|length(unique(temp$Year)) < 2|
       length(unique(temp$SiteYear))<5|length(unique(temp$Reg9Year))<2) {
      mod <- NA
    }else{
      
      if(model == "orig"){
        mod <- try(gam(Total ~ 
                         s(SiteID, bs = "re", k = 5) +
                         s(Reg9Year, bs = "re", k = 5) + 
                         ti(Reg9Year, cumdegday, bs = c("re", "cc"), k = c(5, 10)) +
                         s(cumdegday, bs = "cc", k = 20) +
                         s(Ordinal, bs = "cc", k = 3),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), 
                       # gamma = 1.4, 
                       control = list(maxit = 500)))
        
        
        if("try-error" %in% class(mod)){
          mod <- try(gam(Total ~ 
                           s(SiteID, bs = "re", k = 5) +
                           s(Reg9Year, bs = "re", k = 5) + 
                           # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                           s(cumdegday, bs = "cc", k = 20)+
                           s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
                         control = list(maxit = 500)))
          pars$reduced <- "yes"
        }
      }
      
      
      if(model == "orig.cov"){
        mod <- try(gam(Total ~ 
                         s(zlistlength)+
                         s(ztemperature)+
                         s(zduration)+
                         s(SiteID, bs = "re", k = 5) +
                         s(Reg9Year, bs = "re", k = 5) + 
                         ti(Reg9Year, cumdegday, bs = c("re", "cc"), k = c(5, 10)) +
                         s(cumdegday, bs = "cc", k = 20) +
                         s(Ordinal, bs = "cc", k = 10),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       # method = "REML", 
                       optimizer = c("outer", "newton"), 
                       gamma = 1.4, 
                       control = list(maxit = 500)))
        
        
        if("try-error" %in% class(mod)){
          mod <- try(gam(Total ~ 
                           s(zlistlength)+
                           s(ztemperature)+
                           s(zduration)+
                           s(SiteID, bs = "re", k = 5) +
                           s(Reg9Year, bs = "re", k = 5) + 
                           # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                           s(cumdegday, bs = "cc", k = 20)+
                           s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
                         control = list(maxit = 500)))
          pars$reduced <- "yes"
        }
      }
      
      if(model == "extra"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "cc"), k = c(5, 20), d = c(2, 1)) +
                         s(Ordinal, bs = "cc", k = 10),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), 
                       # gamma = 1.4, 
                       control = list(maxit = 500)))
        
        if("try-error" %in% class(mod)){
          mod <- try(gam(Total ~ 
                           # s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                           s(cumdegday, bs = "cc", k = 20) +
                           s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
                         control = list(maxit = 500)))
          pars$reduced <- "yes"
          
        }
      }
      
      if(model == "extra.cov"){
        mod <- try(gam(Total ~ 
                         s(zlistlength)+
                         s(ztemperature)+
                         s(zduration)+
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "cc"), k = c(5, 20), d = c(2, 1)) +
                         s(Ordinal, bs = "cc", k = 10),
                       # family = nb(theta = NULL, link = "log"),
                       family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), 
                       # gamma = 1.4, 
                       control = list(maxit = 500)))
        
        if("try-error" %in% class(mod)){
          mod <- try(gam(Total ~ 
                           s(zlistlength)+
                           s(ztemperature)+
                           s(zduration)+
                           s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                           s(cumdegday, bs = "cc", k = 20) +
                           s(Ordinal, bs = "cc", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         # gamma = 1.4, 
                         control = list(maxit = 500)))
          pars$reduced <- "yes"
          
        }
      }
      
    }
  }
  
  pars$modtime <- as.period(interval(starttime, Sys.time()))
  pars$AIC <- AIC(mod)
  summod <- summary(mod)
  pars$N <- summod$n
  pars$dev.expl <- summod$dev.expl
  pars$negbin <- mod$family$getTheta(TRUE)
  outlist <- list()
  outlist[["params"]] <- pars
  outlist[["gammod"]] <- mod
  outlist[["datGAM"]] <- temp
  
  #mclust broods
  
  mvmin <- speciesdat$BroodsGAMmin[speciesdat$CommonName==species]
  mvmax <- speciesdat$BroodsGAMmax[speciesdat$CommonName==species]
  
  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday, region9, lat, lon) %>%
    filter(SiteYear %in% unique(temp$SiteYear)) %>% 
    group_by(SiteYear) %>% 
    # filter(yday >= 75 & yday <= 320) %>%
    filter(yday %in% (seq(77, 308, 7) + sample.int(n=6, size=34, replace=TRUE))) %>%
    dplyr::rename(Ordinal = yday)
  pred <- full_join(pred, unique(temp[, c("SiteYear", "Reg9Year")]))
  # pred$zlistlength <- sample(x = temp$zlistlength, size = nrow(pred), replace = TRUE)
  # pred$ztemperature <- sample(x = temp$ztemperature, size = nrow(pred), replace = TRUE)
  # pred$zduration <- sample(x = temp$zduration, size = nrow(pred), replace = TRUE)
  pred$zduration <- 0
  pred$zlistlength <- 0
  pred$ztemperature <- 0
  
  
  # pred$GAM.pred <- as.vector(predict.gam(mod, pred, type = "response", se.fit = TRUE))
  pred <- as.data.frame(pred)
  pred$SiteYear <- as.factor(pred$SiteYear)
  
  outlist[["preddf"]] <- pred
  
  Xp <- predict.gam(object = mod, newdata = pred, type="lpmatrix") ## map coefs to fitted curves
  beta <- coef(mod)
  Vb   <- vcov(mod) ## posterior mean and cov of coefs
  n <- 5 # choose number of simulations
  # set.seed(10)
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  ilink <- family(mod)$linkinv
  linklppreds <- Xp %*% t(mrand)
  nbpreds <- apply(X = linklppreds, 
                   MARGIN = 1, 
                   FUN = function(x){
                     temp <- sort(x)
                     bounds <- quantile(1:n, probs = c(0.025, 0.975))
                     x <- temp[bounds[1]:bounds[2]]
                     x <- ilink(x)
                     x <- rnbinom(n = length(x),
                                  mu = x,
                                  size = mod$family$getTheta(TRUE))
                     return(x)
                   })
  
  
  # #regional
  # SimFunc <- function(nsim){
  simregbroods <- list()
  simsitebroods <- list()
  simpops <- list()
  
  for(nsim in seq_len(nrow(nbpreds))){
    ns <- nsim
    simreglist <- list()
    simsitelist <- list()
    poplist <- list()
    pred$nbpred <- nbpreds[ns,]
    for (j in 1:length(unique(pred$region9))){
      index <- unique(pred$region9)[j]
      clustdat <- pred %>%
        filter(region9 == index)
      
      popdat <- clustdat %>%
        group_by(SiteYear) %>% 
        summarise(totalN = sum(nbpred))
      popdat$nsim <- ns
      poplist[[length(poplist) + 1]] <- popdat
      
      if(sum(clustdat$nbpred) < 10){
        next
      }
      
      dd <- rep(clustdat$cumdegday, clustdat$nbpred)
      daycount <- rep(clustdat$Ordinal, clustdat$nbpred)
      SiteYear <- rep(clustdat$SiteYear, clustdat$nbpred)
      Reg9Year <- rep(clustdat$Reg9Year, clustdat$nbpred)
      dat <- data.frame(SiteYear, Reg9Year, dd, daycount)
      
      mod.dd <- try(Mclust(dat$dd, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
      # mod.day <- try(Mclust(dat$daycount, G=c(mvmin:mvmax), modelNames = "E"), silent = TRUE)
      if("try-error" %in% class(mod.dd) ) next
      # if("try-error" %in% class(mod.day) ) next
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
      simreglist[[length(simreglist) + 1]] <- df1
      
      #attempt SiteYear breakdown of broods
      df2 <- data.frame()
      for (b in 1:length(mod.dd$parameters$mean)){
        brooddat <- dat %>% filter(class==b)
        for (c in 1:length(unique(brooddat$SiteYear))){
          broodtemp <- brooddat %>% filter(SiteYear == unique(brooddat$SiteYear)[c])
          if(length(unique(broodtemp$dd))==1){  #mclust can't fit if only one unique dd value
            df2 <- rbind(df2, data.frame(model = "degday",
                                         nsim = ns,
                                         SiteYear = unique(brooddat$SiteYear)[c],
                                         brood = b,
                                         num = length(broodtemp$dd),
                                         mu = broodtemp$dd[1],
                                         sigma = NA))
          }else{
            broodmod <- try(Mclust(broodtemp$dd, G=1), silent = TRUE)
            if(class(broodmod)=="try-error") next
            df2 <- rbind(df2, data.frame(model = "degday",
                                         nsim = ns,
                                         SiteYear = unique(brooddat$SiteYear)[c],
                                         brood = b,
                                         num = nrow(broodtemp),
                                         mu = broodmod$parameters$mean,
                                         sigma = sqrt(broodmod$parameters$variance$sigmasq)))
          }
        }
      }
      simsitelist[[length(simsitelist) + 1]] <- df2
    }
    
    regoutdf <- data.table::rbindlist(simreglist)
    siteoutdf <- data.table::rbindlist(simsitelist)
    pops <- data.table::rbindlist(poplist)
    
    
    simregbroods[[nsim]] <- regoutdf
    simsitebroods[[nsim]] <- siteoutdf
    simpops[[nsim]] <- pops
    
    # return(outfunc)
  }
  
  
  
  outlist[["simregbroods"]] <- simregbroods
  outlist[["simsitebroods"]] <- simsitebroods
  outlist[["simpops"]] <- simpops
  
  
  return(outlist)
  # saveRDS(outlist, paste(species, model, cutoff, "rds", sep = "."))
}


system.time({
test2 <- mapply(FUN = FitGAM, species = params$species, 
                model = params$model, cutoff = params$cutoff)

})
