rm(list=ls())
source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
# library(pastecs)
library(mclust)
ScaleSumTo1 <- function(x){x/sum(x)}


tpGAM <- function(SiteYearInd, Dat){
  tpdat <- Dat %>%
    filter(SiteYear == SiteYearInd) %>%
    arrange(Ordinal)
  tp <- turnpoints(as.numeric(tpdat$GAM.pred))
  ends <- data.frame(TimeInd = c(1, length(tpdat$GAM.pred)), tptype = c("start", "end"))
  if(length(which(tp$peaks)) >= 1){
    peakpos <- data.frame(TimeInd = which(tp$peaks), tptype = "peak")
    ends <- rbind(ends, peakpos)
  }
  if(length(which(tp$pits)) >= 1){
    pitpos <- data.frame(TimeInd = which(tp$pits), tptype = "pit")
    ends <- rbind(ends, pitpos)
  }
  alltps <- ends %>% arrange(TimeInd)
  alltps$Ordinal <- alltps$TimeInd + min(tpdat$Ordinal) - 1
  alltps$GAM.pred <- tpdat$GAM.pred[alltps$TimeInd]
  alltps$SiteYear <- SiteYearInd
  alltps$aoc <- NA
  alltps$weight <- NA
  for (i in which(alltps$tptype == "peak")){
    alltps$aoc[i] <- sum(tpdat$GAM.pred[alltps$TimeInd[i-1]:alltps$TimeInd[i+1]])
    alltps$weight[i] <- alltps$aoc[i] / sum(tpdat$GAM.pred)
  }
  return(alltps)
}


# try with any species

data <- fread("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/data.trim.csv", header = TRUE)
# data <- data.table(data)
setnames(data,"SiteID.x","SiteID")
data[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]
data$SiteDate <- lubridate::ymd(as.character(data$SiteDate))
data$CommonName <- as.character(data$CommonName)
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week"), with = FALSE])
# data <- data[, list(SeqID, SiteID.x, SiteDate, Week, Total, CheckListKey, CommonName)]


site_geo <- read.csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled_update2016.csv", 
                     header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

####
# cluster sites into regions, so that RegionYear can be factor used in GAMs
# this reduces # of random effects and speeds up model fitting

# sites <- read.csv("data/OHsites_reconciled.csv")
sites <- read.csv("data/OHsites_reconciled_update2016.csv")

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

sites$region4 <- mapvalues(sites$region4, from = c("1", "2", "3", "4"), 
                           to = c("NE", "NW", "CN", "SW"))
sites$region9 <- mapvalues(sites$region9, from = c(as.character(1:9)),
                           to = c("Cuyahoga", "Painesville", "Huron",
                                  "Findlay", "Columbus", "Cincinnati",
                                  "Toledo", "Dayton", "Lancaster"))
sites$SiteID <- formatC(sites$site, width = 3, format = "d", flag = "0")

# plot(sitemod, what = "density")

# visualize clusters
sites$class <- as.character(sitemod$classification)
a <- ggplot(data = sites, aes(x = lon, y = lat, group = class, color = class)) + geom_point()
a
site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))




allspec <- unique(data$CommonName)
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspec[1:122])),
            temperature = mean(c(StartTemp[1], EndTemp[1]), na.rm = TRUE),
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


gdd <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/daymet/growingDD_Daymet.RDS")
gdd <- gdd[gdd$year >= 1995, ]
gdd <- as.data.frame(gdd)
gdd$date <- strptime(paste(gdd$year, gdd$yday, sep = " "), format = "%Y %j")
gdd$SiteDate <- as.Date(gdd$date)
gdd$date <- NULL

gdd <- gdd %>%
  mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))

SpeciesList <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/SpeciesList.rds")

species <- SpeciesList %>% arrange(Present) %>% select(CommonName)


# gamlist <- as.list(species$CommonName)
# for (i in 1:nrow(species)){
  GAMspecies <- function(i){
  sp <- species[i, ]
  gamlist <- list(sp)
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
  counts$listlength <- scale(counts$listlength)
  counts$temperature <- scale(counts$temperature)
  counts$duration <- scale(counts$duration)
  counts$SiteDate <- as.Date(counts$SiteDate)
  
  print(sp)
  
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  
  datGAM <- counts[YearTotal >= 3]
  datGAM <- datGAM[SurvPerYear >= 15]
  if(nrow(datGAM) == 0) next
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
      filter(YearTotal > 1,
             SurvSeen > 0) %>%
      dplyr::select(SiteID, Total, listlength, Ordinal, Year, temperature, duration) %>% 
      as.data.frame()
    
    # pad zeros at beginning and end for GAM fitting
    zeros <- c(50, 60, 330, 340)
    tempdf <- expand.grid(unique(temp$SiteID), zeros)
    names(tempdf) <- c("SiteID", "Ordinal")
    tempdf$SiteID <- as.character(tempdf$SiteID)
    tempdf$listlength <- 0
    tempdf$temperature <- 0
    tempdf$duration <- 0
    tempdf$Total <- 0

    test <- merge(tempdf, distinct(dplyr::select(temp, SiteID, Year)),
                  by = "SiteID", all.x = TRUE, all.y = FALSE)
    
    outlist[[y-min(unique(dat$Year))+1]] <- rbind(temp, test)
  }
  
  dat <- rbindlist(outlist)
  dat$Year <- as.factor(as.character(dat$Year))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
  dat$temperature[which(is.na(dat$temperature))] <- 0
  dat$duration[which(is.na(dat$duration))] <- 0
  dat <- merge(dat, site_geo, by = c("SiteID"))
  dat$Reg9Year <- as.factor(paste(dat$region9, dat$Year, sep = "_"))
  dat$Reg4Year <- as.factor(paste(dat$region4, dat$Year, sep = "_"))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$region9 <- as.factor(dat$region9)
  
  dat <- dat %>%
    group_by(region9) %>%
    mutate(uniqSiteYear = length(unique(SiteYear))) %>%
    filter(uniqSiteYear > 1)
  
  dat$SiteDate <- as.Date(paste(dat$Ordinal, dat$Year, sep = "-"),
                          format = "%j-%Y")
  # trying to add GDD instead of ordinal date
  dat <- merge(dat, gdd, by = c("SiteID", "SiteDate"), all.x = TRUE, all.y = FALSE)
  dat$SiteDate <- as.character(dat$SiteDate)
  
  temp <- dat
  temp <- temp %>% 
    # filter(Year %in% c(2004, 2005, 2006, 2007)) %>% 
    dplyr::select(SiteID, Total, Ordinal, Year, SiteYear, lat, lon,
           region9, region4, Reg9Year, Reg4Year, cumdegday) %>% 
    data.frame()
  temp$SiteID <- as.factor(temp$SiteID)
  temp <- droplevels.data.frame(temp)
  
  
  if(sum(temp$Total) < 20) next
  if(length(unique(temp$SiteID)) < 2) next
  if(length(unique(temp$Year)) < 2) next
  

  
  # again, got cold feet and went back and tried lots of gam formula
  # mod7 was what i used for mvmod, and it's the best by AIC (for NPE)
  # stop playing around with this!!!
  # system.time({mod5 <- try(mod5 <- gam(Total ~ 
  #                                        s(SiteID, bs = "re") +
  #                                        region9 +
  #                                          s(cumdegday, k = 20, by = region9) +
  #                                          s(cumdegday, Year, bs = "fs", m=1),
  #                                        family = nb(theta = NULL, link = "log"),
  #                                        # family = poisson(link = "log"),
  #                                        data = temp,
  #                                      method = "REML", 
  #                                      optimizer = c("outer", "newton"), 
  #                                        gamma = 1.4, 
  #                                        control = list(maxit = 500)))
  # })
  # 
  # system.time({mod5b <- try(mod5b <- gam(Total ~ 
  #                                        SiteID + 
  #                                        Reg9Year +
  #                                        s(cumdegday, k = 20, by = Reg9Year),
  #                                      family = nb(theta = NULL, link = "log"),
  #                                      # family = poisson(link = "log"),
  #                                      data = temp,
  #                                      method = "GCV.Cp", 
  #                                      optimizer = "perf", 
  #                                      gamma = 1.4, 
  #                                      control = list(maxit = 500)))
  # })
  # 
  # system.time({mod6 <- try(mod6 <- gam(Total ~ 
  #                                          SiteID +
  #                                          Year +
  #                                          s(cumdegday, Reg9Year, bs = "fs", m=1),
  #                                        family = nb(theta = NULL, link = "log"),
  #                                        # family = poisson(link = "log"),
  #                                        data = temp,
  #                                      method = "REML", 
  #                                      optimizer = c("outer", "newton"), 
  #                                        gamma = 1.4, 
  #                                        control = list(maxit = 500)))
  # })
  # system.time({mod6b <- try(mod6b <- gam(Total ~ 
  #                                        s(SiteID, bs = "re") +
  #                                        s(Year, bs = "re") +
  #                                        s(cumdegday, Reg9Year, bs = "fs", m=1),
  #                                      family = nb(theta = NULL, link = "log"),
  #                                      # family = poisson(link = "log"),
  #                                      data = temp,
  #                                      method = "REML", 
  #                                      optimizer = c("outer", "newton"), 
  #                                      gamma = 1.4, 
  #                                      control = list(maxit = 500)))
  # })
  
  # system.time({mod7 <- try(mod7 <- gam(Total ~ 
  #                           s(SiteID, bs = "re", k = 5) +
  #                           s(Reg9Year, bs = "re", k = 5) + 
  #                           ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
  #                           s(cumdegday, bs = "cr", k = 20),
  #                         family = nb(theta = NULL, link = "log"),
  #                         # family = poisson(link = "log"),
  #                         data = temp,
  #                         method = "REML", 
  #                         optimizer = c("outer", "newton"), 
  #                         gamma = 1.4, 
  #                         control = list(maxit = 500)))
  # })
  
  system.time({mod7b <- try(mod7b <- gam(Total ~ 
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
  
  # system.time({mod8 <- try(gam(Total ~ s(listlength, k = 3) +  
  #                   s(temperature, k = 3) +
  #                   s(duration, k = 3) +
  #                   s(SiteID, bs = "re") + 
  #                   s(Year, bs = "re") +
  #                   s(Year, cumdegday, bs = "fs", m = 1) +
  #                   te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
  #                 family = nb(theta = NULL, link = "log"),
  #                 # family = poisson(link = "log"),
  #                 data = temp,
  #                 method = "GCV.Cp", 
  #                 optimizer = "perf", 
  #                 gamma = 1.4, 
  #                 control = list(maxit = 500)))
  # })
#   system.time({
#   mod9 <- try(gam(Total ~ #s(listlength, k = 3) +  
#                     #s(temperature, k = 3) +
#                     #s(duration, k = 3) +
#                     s(SiteID, bs = "re") + 
#                     te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1), by = Year),
#                   family = nb(theta = NULL, link = "log"),
#                   # family = poisson(link = "log"),
#                   data = temp,
#                   method = "REML", 
#                   optimizer = c("outer", "newton"), 
#                   gamma = 1.4, 
#                   control = list(maxit = 500)))
# })
  # 
  # if(("gam" %in% class(mod7)) == FALSE){
  #   mod7 <- try(mod7 <- gam(Total ~ 
  #                             s(SiteID, bs = "re", k = 5) +
  #                             s(Reg9Year, bs = "re", k = 5) + 
  #                             # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
  #                             s(cumdegday, bs = "cr", k = 20),
  #                           family = nb(theta = NULL, link = "log"),
  #                           # family = poisson(link = "log"),
  #                           data = temp,
  #                           method = "REML", 
  #                           optimizer = c("outer", "newton"), 
  #                           sp = c(-1, -1, -1),
  #                           gamma = 1.4, 
  #                           control = list(maxit = 500)))
  # }
  
  if(("gam" %in% class(mod7b)) == FALSE){
    mod7b <- try(mod7b <- gam(Total ~ 
                              s(SiteID, bs = "re", k = 5) +
                              s(Reg9Year, bs = "re", k = 5) + 
                              # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                              s(cumdegday, bs = "cr", k = 20) +
                                s(Ordinal, bs = "cr", k = 10),
                            family = nb(theta = NULL, link = "log"),
                            # family = poisson(link = "log"),
                            data = temp,
                            method = "REML", 
                            optimizer = c("outer", "newton"), 
                            sp = c(-1, -1, -1),
                            gamma = 1.4, 
                            control = list(maxit = 500)))
  }
  
  
  if(("gam" %in% class(mod7b)) == FALSE) next
  
  # 
  # temp1 <- temp %>% filter(Ordinal > 70 & Ordinal < 320)
  # # 
  # # start.surv <- sort(unique(temp1$Ordinal))[1]
  # # end.surv <- sort(unique(temp1$Ordinal), decreasing = TRUE)[1]
  # # allDays <- expand.grid(c(start.surv:end.surv), unique(temp1$SiteYear))
  # # names(allDays) <- c("Ordinal", "SiteYear")
  # # t <- str_split_fixed(allDays$SiteYear, pattern = "_", 2)
  # # allDays$SiteID <- t[,1]
  # # allDays$Year <- t[,2]
  # 
  # start.gdd <- sort(unique(temp1$cumdegday))[1]
  # end.gdd <- sort(unique(temp1$cumdegday), decreasing = TRUE)[1]
  # allgdd <- expand.grid(seq(start.gdd, end.gdd, by = 25), unique(temp1$SiteYear))
  # names(allgdd) <- c("cumdegday", "SiteYear")
  # t <- str_split_fixed(allgdd$SiteYear, pattern = "_", 2)
  # allgdd$SiteID <- t[,1]
  # allgdd$Year <- t[,2]
  # 
  # 
  # newData <- allgdd
  # # newData <- allDays
  # geo <- temp1 %>% dplyr::select(SiteID, Ordinal, lat, lon, Reg9Year, region9, Year) %>% distinct()
  # newData <- merge(newData, geo, by = c("SiteID", "Year"), all.x = TRUE, all.y = FALSE)
  # # newData$listlength <- mean(temp1$listlength, na.rm = TRUE)
  # # newData$temperature <- mean(temp1$temperature, na.rm = TRUE)
  # # newData$duration <- mean(temp1$duration, na.rm = TRUE)
  # # 
  # temp$GAM.pred <- as.vector(predict.gam(mod7, temp, type = "response"))
  # temp$GAM.pred.b <- as.vector(predict.gam(mod7b, temp, type = "response"))
  # 
  # # temp$GAM.pred <- predict.gam(mod7, type = "response")
  # # newData$GAM.pred <- as.vector(predict.gam(mod7, newData, type = "response",
  # #                                 exclude = "ti(Reg9Year,cumdegday"))
  # newData <- data.table(temp)
  # newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"] 
  # newData[, Gamma.b := ScaleSumTo1(GAM.pred.b), by = "SiteYear"] 
  # 
  # newData$SiteID <- as.factor(newData$SiteID)
  # newData$Year <- as.factor(newData$Year)
  # newData <- as.data.frame(newData)
  # # newData <- newData %>%
  # #   group_by(SiteYear) %>% 
  # #   mutate(sumEst = sum(GAM.pred)) %>% 
  # #   filter(sumEst >= 20)
  #   
  # 
  # c <- ggplot(data = newData, aes(x = cumdegday, y = Gamma.b, group = SiteID, color = SiteID)) +
  #   geom_line(size = 1, alpha = .5) + 
  #   theme_bw() + theme(legend.position = "none") +
  #   facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
  # print(c)
  # d <- ggplot(data = newData, aes(x = cumdegday, y = GAM.pred.b, group = SiteID, color = SiteID)) + 
  #   geom_line(size = 1, alpha = .5) +
  #   theme_bw() + geom_point(data = temp1, aes(x = cumdegday, y = Total))
  # d <- d + facet_wrap( ~ Year, scales = "free_y") + theme(legend.position = "none") +
  #   ggtitle(paste("GAM Predictions", sp, sep = " "))
  # print(d)
  # 
  # 
  # pdf(paste("CovGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
  # print(c)
  # dev.off()
  # 
  # pdf(paste("CovGAMPRED", sp, ".pdf", sep = ""), width = 10, height = 6)
  # print(d)
  # dev.off()
  
  # gamlist$mod <- mod7
  gamlist$modb <- mod7b
  gamlist$counts <- temp

  saveRDS(gamlist, paste("gam.gdd.plus.day", sp, ".rds", sep = ""))
}


i <- 53
sp <- species[i, "CommonName"]
gamlist <- readRDS(paste("gam", sp, ".rds", sep = ""))


# errors: running out of memory???
 
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
                          "species"))
  
  # test <- parLapply(cl, c(105, 107, 108, 109, 112, 114, 116), GAMspecies)
  test <- parLapply(cl, c(104, 106, 110, 111, 113, 115, 117,118,119,120,121,122), GAMspecies)
  
   stopCluster(cl)
})


test <- lapply(as.list(c(121,122)),
               GAMspecies)

saveRDS(test, file = "BroodMixSigLim2.rds")


