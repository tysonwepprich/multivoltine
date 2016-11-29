
#final!
#1-80 species at 4 cpus per node
#461219, slr 7580
#81-100 species at 2 cpus per node
#461225, slr 7735
#ALL WRAPPER ERRORS!?!

#then check for wrapper code error, possibly memory related if too parallel
#from slr7580, reran these as
#461231, slr7035
#worked with no wrapper code errors

#trying again after errors with starttime. Trying first 100 species

source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
library(mclust)
ScaleSumTo1 <- function(x){x/sum(x)}



# try with any species

data <- fread("data/data.trim.csv", header = TRUE)
# data <- data.table(data)
setnames(data,"SiteID.x","SiteID")
data[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]
data$SiteDate <- lubridate::ymd(as.character(data$SiteDate))
data$CommonName <- as.character(data$CommonName)
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week"), with = FALSE])
# data <- data[, list(SeqID, SiteID.x, SiteDate, Week, Total, CheckListKey, CommonName)]


site_geo <- fread("data/OHsites_reconciled_update2016.csv", header = TRUE)
# site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

allspec <- unique(data$CommonName)
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspec[1:122])),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct()

#####
# Aside, PCA of covariates
# covdata$temperature[which(covdata$temperature < 50)] <- NA
# covdata$duration[which(covdata$duration == 0)] <- NA
# covdata$temperature <- scale(covdata$temperature)
# covdata$listlength <- scale(covdata$listlength)
# covdata$duration <- scale(covdata$duration)
# covdata$temperature[which(is.na(covdata$temperature))] <- 0
# covdata$duration[which(is.na(covdata$duration))] <- 0
# 
# covpca <- prcomp(covdata[,c(2:4), with = FALSE])

#####

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

# visualize clusters
# sites$class <- as.character(sitemod$classification)
# a <- ggplot(data = sites, aes(x = lon, y = lat, group = region9, color = region9)) + geom_point()
# a

site_geo <- fread("data/OHsites_reconciled_update2016.csv")
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
gdd_tomerge <- gdd %>% dplyr::select(SiteYear, yday, cumdegday) %>%
  dplyr::rename(Ordinal = yday)




# gdd <- readRDS("data/growingDD_Daymet.RDS")
# gdd <- gdd[gdd$year >= 1995, ]
# gdd <- as.data.frame(gdd)
# gdd$date <- strptime(paste(gdd$year, gdd$yday, sep = " "), format = "%Y %j")
# gdd$SiteDate <- as.Date(gdd$date)
# 
# gdd <- gdd %>%
#   mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))

SpeciesList <- readRDS("data/SpeciesList.rds")

species <- SpeciesList %>% arrange(Present) %>% select(CommonName)

# models <-  c("mod1", "mod1simp", "mod2", "mod2simp", "mod1cov", "mod1covsimp", "mod2cov", "mod2covsimp", "overparam")
models <- c("orig", "extra")
cutoff <- c("strict", "loose")
params <- expand.grid(species$CommonName, models, cutoff,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model", "cutoff")

#testing
# params <- params %>% filter(species == "Carolina Satyr")
#redoing wrapper errors
ok <- redo %>% select(species, model, cutoff)
params <- anti_join(params, ok)
params <- params[c(17:20),]


params <- params[sample(1:nrow(params)), ] #rearrange for parallel comp speed

# # data_file Rdata
# dataIN <- c("gdd", "data", "surveys", "covdata", "site_geo", "params")
# save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
# paramIN <- data.frame(nRow = seq(1:nrow(params)))

# test <- apply(X = params, MARGIN = 1, FitGAM)

# calculate null hypotheses for M for different species
# gams <- slurm_apply(f = FitGAM, params = paramIN, 
#                          cpus_per_node = 8, nodes = 4, 
#                          data_file = "dataIN.RData", 
#                          # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
#                          output = "raw")
ty <- slurm_apply(f = FitGAM, params = params, nodes = 1,
                  cpus_per_node = 4, 
                  add_objects = c("gdd", "data", "surveys", "covdata", "site_geo", "params"),
                  slurm_options = list(partition = "sesynctest"))

test <- mapply(FUN = FitGAM, species = params$species, model = params$model, cutoff = params$cutoff)



# gamlist <- as.list(species$CommonName)
for (i in 1:nrow(params)){
  sp <- species <- params$species[i]
  model <- params$model[i]
  cutoff <- params$cutoff[i]
# FitGAM <- function(nRow){
# FitGAM <- function(species, model, cutoff){
  # sp <- species
  # model <- model
  # cutoff <- cutoff
  reduced <- NA
  pars <- data.frame(species, model, cutoff, reduced)
  # pars <- params[nRow,]
  # sp <- pars$species
  # model <- pars$model
  counts <- data %>% filter(CommonName == sp) %>% data.table()
  counts[, Ordinal := yday(SiteDate)]
  counts[, Year := year(SiteDate)]
  
  
  survs <- surveys[year(SiteDate) %in% unique(counts$Year)]
  survs <- survs[SiteID %in% unique(counts$SiteID)]
  survs$Year <- year(survs$SiteDate)
  #Add zeros to surveys when species not counted during a survey
  
  test <- merge(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"), all.x = TRUE)
  counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year"), with = FALSE]
  counts$Total <- mapvalues(counts[,Total], from = NA, to = 0)
  counts <- merge(counts, site_geo, by = "SiteID")
  counts <- merge(counts, covdata, by = "SeqID", all.x = TRUE, all.y = FALSE)
  counts$temperature[which(counts$temperature < 50)] <- NA
  counts$duration[which(counts$duration == 0)] <- NA
  
  # #scaled over whole season/state
  # counts$listlength <- scale(counts$listlength)
  # counts$temperature <- scale(counts$temperature)
  # counts$duration <- scale(counts$duration)
  #scaled by region9/week
  counts <- counts[, `:=` (zlistlength = as.numeric(scale(listlength)),
                           ztemperature = as.numeric(scale(temperature)),
                           zduration = as.numeric(scale(duration))),
                   by = list(region9, Week)]
                            
  
  # trying to add GDD instead of ordinal date
  counts <- merge(counts, gdd, by = c("SiteID", "SiteDate", "lat", "lon", "region9", "region4"), all.x = TRUE, all.y = FALSE)
  
  print(sp)
  
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
  starttime <- Sys.time()
  
  if(nrow(datGAM) == 0){
    mod <- NA
  }else{
    datGAM[, SitesObserved := length(unique(SiteID)), by = list(Year)]
    datGAM <- datGAM[SitesObserved >= 1]
    
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
        dplyr::select(SiteID, Total, listlength, lat, lon, region9, Ordinal, Year, temperature, duration, cumdegday) %>% 
        data.frame()
      
      # pad zeros at beginning and end for GAM fitting
      zeros <- c(50, 60, 330, 340)
      tempdf <- expand.grid(unique(temp$SiteID), zeros)
      names(tempdf) <- c("SiteID", "yday")
      tempdf$year <- y
      tempdf$listlength <- 0
      tempdf$temperature <- 0
      tempdf$duration <- 0
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
    dat$temperature[which(is.na(dat$temperature))] <- 0
    dat$duration[which(is.na(dat$duration))] <- 0
    dat$Reg9Year <- as.factor(paste(dat$region9, dat$Year, sep = "_"))
###
    #this chunk was in sesync code, but not others for fitgam!
    #throws out some SiteYears I didn't expect
    
    dat <- dat %>%
      group_by(region9) %>%
      mutate(uniqSiteYear = length(unique(SiteYear))) %>%
      filter(uniqSiteYear > 1)
    ###
    
    
    temp <- dat

    if(sum(temp$Total) < 20|length(unique(temp$SiteID)) < 2|length(unique(temp$Year)) < 2|
       length(unique(temp$SiteYear))<5|length(unique(temp$Reg9Year))<2) {
      mod <- NA
    }else{
      

      if(model == "orig"){
        mod <- try(gam(Total ~ 
                         s(SiteID, bs = "re", k = 5) +
                         s(Reg9Year, bs = "re", k = 5) + 
                         ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                         s(cumdegday, bs = "cr", k = 20) +
                         s(Ordinal, bs = "cr", k = 10),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), 
                       gamma = 1.4, 
                       control = list(maxit = 500)))

        
        if(class(mod)=="try-error"){
          mod <- try(gam(Total ~ 
                           s(SiteID, bs = "re", k = 5) +
                           s(Reg9Year, bs = "re", k = 5) + 
                           # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                           s(cumdegday, bs = "cr", k = 20)+
                           s(Ordinal, bs = "cr", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), 
                         gamma = 1.4, 
                         control = list(maxit = 500)))
          pars$reduced <- "yes"
        }
      }
      
      
      
      if(model == "extra"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 20), d = c(2, 1)) +
                         s(Ordinal, bs = "cr", k = 10),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))

        if(class(mod)=="try-error"){
          mod <- try(gam(Total ~ 
                           s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                           s(cumdegday, bs = "tp", k = 20) +
                           s(Ordinal, bs = "cr", k = 10),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))
          pars$reduced <- "yes"
          
          }
      }
      
      if(model == "region"){
        temp$region9 <- as.factor(temp$region9)
        mod <- try(gam(Total ~ 
                         s(SiteYear, bs = "re", k = 5) +
                         s(cumdegday, bs = "cr", k = 20, by = region9),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), 
                       gamma = 1.4, 
                       control = list(maxit = 500)))
        if(class(mod)=="try-error"){
          mod <- try(gam(Total ~ 
                           s(SiteYear, bs = "re", k = 5) +
                           s(cumdegday, bs = "cr", k = 20),
                         family = nb(theta = NULL, link = "log"),
                         # family = poisson(link = "log"),
                         data = temp,
                         method = "REML", 
                         optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))
          pars$reduced <- "yes"
          
        }
      }
      
      
    } 
  }
  
  pars$modtime <- as.numeric(Sys.time() - starttime)
  outlist <- list()
  outlist[[1]] <- pars
  outlist[[2]] <- mod
  # return(outlist)
  saveRDS(outlist, paste(species, model, cutoff, "rds", sep = "."))
}


library(parallel)
# multicore
system.time({
  cl <- makeCluster(2)
  clusterEvalQ(cl, {
    library(mclust)
    library(data.table)
    library(plyr)
    library(dplyr)
    library(mgcv)
  })
  clusterExport(cl=cl,
                varlist=c("site_geo",
                          "gdd",
                          "surveys",
                          "data",
                          "covdata",
                          "params"))
  
  # test <- parLapply(cl, c(105, 107, 108, 109, 112, 114, 116), GAMspecies)
  # test <- parLapply(cl, c(104, 106, 110, 111, 113, 115, 117,118,119,120,121,122), FitGAM)
  # params <- as.array(params, dimnames=list("species", "model"))
  # dimnames(params) <- list("species", "model")
  test <- parApply(cl, params, 1, FitGAM)
  
  stopCluster(cl)
})

system.time({test <- mapply(FUN = FitGAM, species = params$species, model = params$model,
                            cutoff = params$cutoff)
})


# get results
dat <- readRDS("_rslurm_slr8547/results_1.RData")


# extract data from SlurmCov results
slurm_codes <- c("_rslurm_slr7735")
slurm_out <- list()
outlist <- list()
# setwd("slurmCovOutput")

for (j in 1:length(slurm_codes)){
  missing_files <- c()
  tmpEnv <- new.env()
  for (i in 0:12) {
    fname <- paste0("results_", i, 
                    ".RData")
    if(fname %in% list.files(slurm_codes[j])){
      slurm_out <- readRDS(paste0(slurm_codes[j], "/", fname))
    }else{
      next
    }
    for (k in 1:length(slurm_out)){
      outdf <- slurm_out[[k]][[1]]
      if(class(outdf) == "character"){
        outdf <- data.frame(species = names(slurm_out)[k],
                            model = NA,
                            cutoff = NA,
                            modtime = NA,
                            AIC = NA,
                            error = "wrapper code",
                            cluster = i,
                            njob = k)
      }else{
        if(is.na(slurm_out[[k]][[2]])){
          outdf <- data.frame(species = names(slurm_out)[k],
                              model = slurm_out[[k]][[1]]$model,
                              cutoff = slurm_out[[k]][[1]]$cutoff,
                              modtime = NA,
                              AIC = NA,
                              error = "not enough data",
                              cluster = i,
                              njob = k)
        }else{
          if(class(slurm_out[[k]][[2]]) == "try-error"){
            outdf <- data.frame(species = names(slurm_out)[k],
                                model = NA,
                                cutoff = NA,
                                modtime = NA,
                                AIC = NA,
                                error = "gam try-error",
                                cluster = i,
                                njob = k)
          }else{
            outdf$AIC <- AIC(slurm_out[[k]][[2]])
            outdf$error <- "ok"
            outdf$cluster <- i
            outdf$njob <- k
          }
        }
      }
      outlist[[length(outlist)+1]] <- outdf
    }
  }
}

test <- bind_rows(outlist)

test %>% arrange(cluster, njob ,species, cutoff, model) %>% data.frame()

test2 <- test %>% 
  group_by(error, cluster, species, model) %>% data.frame()

redo <- test %>% filter(error != "wrapper code")


#graveyard


system.time({mod <- try(gam(Total ~ 
                 s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
                 te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
               family = nb(theta = NULL, link = "log"),
               # family = poisson(link = "log"),
               data = temp,
               method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
})

system.time({mod2 <- try(mod7b <- gam(Total ~ 
                                         s(SiteID, bs = "re", k = 5) +
                                         s(Reg9Year, bs = "re", k = 5) + 
                                         ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                                         s(cumdegday, bs = "cr", k = 20) +
                                         s(Ordinal, bs = "cr", k = 10),
                                       family = nb(theta = NULL, link = "log"),
                                       # family = poisson(link = "log"),
                                       data = temp,
                                       method = "REML", 
                                       optimizer = c("outer", "newton"), 
                                       gamma = 1.4, 
                                       control = list(maxit = 500)))
})


system.time({mod3 <- try(mod7b <- gam(Total ~ 
                                         s(SiteID, bs = "re", k = 5) +
                                         s(Year, bs = "re", k = 5) + 
                                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)) +
                                         s(Ordinal, bs = "cr", k = 10),
                                       family = nb(theta = NULL, link = "log"),
                                       # family = poisson(link = "log"),
                                       data = temp,
                                       method = "REML", 
                                       optimizer = c("outer", "newton"), 
                                       gamma = 1.4, 
                                       control = list(maxit = 500)))
})


system.time({mod4 <- try(gam(Total ~ s(SiteID, bs = "re") + s(Year, bs = "re") +
                              s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
                              te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                            family = nb(theta = NULL, link = "log"),
                            # family = poisson(link = "log"),
                            data = temp,
                            method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
})


system.time({mod5 <- try(gam(Total ~ s(SiteID, bs = "re") + s(Year, bs = "re") +
                               s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                               te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)) +
                               s(Ordinal, bs = "cr", k = 10),
                             family = nb(theta = NULL, link = "log"),
                             # family = poisson(link = "log"),
                             data = temp,
                             method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
})

AIC(mod, mod2, mod3, mod4, mod5)


gammod <- test[[2]]
gammod <- slurm_out[[1]][[2]]

datGAM <- temp
pred <- gdd %>%
  dplyr::select(SiteID, SiteYear, yday, cumdegday, lat, lon, region9, year, meanGDD) %>%
  filter(SiteYear %in% unique(datGAM$SiteYear)) %>% 
  filter(yday >= 50 & yday <= 340) %>% 
  dplyr::rename(Ordinal = yday,
                Year = year)
pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year"), with = FALSE]))

pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
# pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response",
#   exclude = c("s(SiteID)")  ))
pred <- pred %>% 
  group_by(SiteYear) %>% 
  mutate(SiteYearGDD = max(cumdegday),
         Gamma = GAM.pred.reg9yr / sum(GAM.pred.reg9yr))
pred <- as.data.frame(pred)

c <- ggplot(data = pred, aes(x = cumdegday, y = Gamma, group = SiteID, color = SiteID)) +
  geom_line(size = 1, alpha = .5) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
print(c)

d <- ggplot(data = temp, aes(x = cumdegday, y = Total, group = SiteID, color = SiteID)) +
  geom_point(size = 1, alpha = .5) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
print(d)

