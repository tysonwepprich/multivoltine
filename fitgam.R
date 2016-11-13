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
params <- expand.grid(species$CommonName, models,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "model")

params <- params[c(70:73, 192:195),]
# params <- params[sample(1:nrow(params)), ] #rearrange for parallel comp speed

# data_file Rdata
dataIN <- c("gdd", "data", "surveys", "covdata", "site_geo", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRow = seq(1:nrow(params)))

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

test <- mapply(FUN = FitGAM, species = params$species, model = params$model)
# gamlist <- as.list(species$CommonName)
# for (i in 1:nrow(species)){
# FitGAM <- function(nRow){
FitGAM <- function(species, model){
  sp <- species
  model <- model
  pars <- data.frame(species, model)
  # pars <- params[nRow,]
  # sp <- pars$species
  # model <- pars$model
  counts <- data[CommonName == sp]
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
  counts$listlength <- scale(counts$listlength)
  counts$temperature <- scale(counts$temperature)
  counts$duration <- scale(counts$duration)
  
  # trying to add GDD instead of ordinal date
  counts <- merge(counts, gdd, by = c("SiteID", "SiteDate", "lat", "lon", "region9", "region4"), all.x = TRUE, all.y = FALSE)
  
  print(sp)
  
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  
  datGAM <- counts[YearTotal >= 1]
  datGAM <- datGAM[SurvPerYear >= 10]
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
      names(tempdf) <- c("SiteID", "Ordinal")
      tempdf$listlength <- 0
      tempdf$temperature <- 0
      tempdf$duration <- 0
      tempdf$Total <- 0
      tempdf$cumdegday <- c(rep(0, nrow(tempdf)/2), rep(max(temp$cumdegday), nrow(tempdf)/2))
      
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
    
    temp <- dat
    
    if(sum(temp$Total) < 20|length(unique(temp$SiteID)) < 2|length(unique(temp$Year)) < 2) {
      mod <- NA
    }else{
      
      
      if(model == "mod1"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "mod1simp"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) + 
                         s(Ordinal, bs = "tp", k = 15),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp, 
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      
      if(model == "mod2"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if (model == "mod2simp"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) + 
                         s(cumdegday, bs = "tp", k = 15),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp, 
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "mod1cov"){
        mod <- try(gam(Total ~ s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                         s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "mod1covsimp"){
        mod <- try(gam(Total ~ s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                         s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) + 
                         s(Ordinal, bs = "tp", k = 15),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp, 
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "mod2cov"){
        mod <- try(gam(Total ~ s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "mod2covsimp"){
        mod <- try(gam(Total ~ s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) + 
                         s(cumdegday, bs = "tp", k = 15),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp, 
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
      if(model == "overparam"){
        mod <- try(gam(Total ~ 
                         s(SiteYear, cumdegday, bs = "fs", k = 10, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(7, 20), d = c(2, 1)),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
      }
      
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
      }
      
      if(model == "extra"){
        mod <- try(gam(Total ~ s(SiteID, bs = "re") + s(Year, bs = "re") +
                         s(SiteYear, cumdegday, bs = "fs", k = 5, m = 1) +
                         te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)) +
                         s(Ordinal, bs = "cr", k = 10),
                       family = nb(theta = NULL, link = "log"),
                       # family = poisson(link = "log"),
                       data = temp,
                       method = "REML", 
                       optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))
      }
    } 
  }
  
  outlist <- list()
  outlist[[1]] <- pars
  outlist[[2]] <- mod
  return(outlist)
}

library(parallel)
# multicore
system.time({
  cl <- makeCluster(4)
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