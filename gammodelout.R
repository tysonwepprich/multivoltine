source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
library(pastecs)
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


site_geo <- read.csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled.csv", header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

####
# cluster sites into regions, so that RegionYear can be factor used in GAMs
# this reduces # of random effects and speeds up model fitting

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

sites$region4 <- mapvalues(sites$region4, from = c("1", "2", "3", "4"), 
                           to = c("NE", "NW", "CN", "SW"))
sites$region9 <- mapvalues(sites$region9, from = c(as.character(1:9)),
                           to = c("Cuyahoga", "Painesville", "Huron",
                                  "Findlay", "Columbus", "Cincinnati",
                                  "Toledo", "Dayton", "Lancaster"))
sites$SiteID <- formatC(sites$site, width = 3, format = "d", flag = "0")

# plot(sitemod, what = "density")

# visualize clusters
# sites$class <- as.character(sitemod$classification)
# a <- ggplot(data = sites, aes(x = lon, y = lat, group = class, color = class)) + geom_point()
# a
site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))




allspec <- unique(data$CommonName)
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspec[1:122])),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration) %>%
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

gdd <- gdd %>%
  mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))

SpeciesList <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/SpeciesList.rds")

species <- SpeciesList %>% arrange(Present) %>% select(CommonName)


# gamlist <- as.list(species$CommonName)
for (i in 116:nrow(species)){
  sp <- species[i, ]
  gamlist <- list(sp)
  counts <- data %>% filter(CommonName == sp)
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
  
  
  # trying to add GDD instead of ordinal date
  counts <- merge(counts, gdd, by = c("SiteID", "SiteDate"), all.x = TRUE, all.y = FALSE)
  
  print(sp)
  
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  
  datGAM <- counts[YearTotal > 1]
  datGAM <- datGAM[SurvPerYear >= 10]
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
      select(SiteID, Total, listlength, lat, lon, Ordinal, Year, temperature, duration, cumdegday)
    
    # pad zeros at beginning and end for GAM fitting
    zeros <- c(50, 60, 330, 340)
    tempdf <- expand.grid(unique(temp$SiteID), zeros)
    names(tempdf) <- c("SiteID", "Ordinal")
    tempdf$listlength <- 0
    tempdf$temperature <- 0
    tempdf$duration <- 0
    tempdf$Total <- 0
    tempdf$cumdegday <- c(rep(0, nrow(tempdf)/2), rep((max(temp$cumdegday) + 200), nrow(tempdf)/2))
    
    test <- merge(tempdf, distinct(select(temp, SiteID, lat, lon, Year)),
                  by = "SiteID", all.x = TRUE, all.y = FALSE)
    
    outlist[[y]] <- rbind(temp, test)
  }
  
  dat <- rbindlist(outlist)
  dat$Year <- as.factor(as.character(dat$Year))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
  dat$temperature[which(is.na(dat$temperature))] <- 0
  dat$duration[which(is.na(dat$duration))] <- 0
  dat <- merge(dat, site_geo, by = c("SiteID", "lat", "lon"))
  dat$Reg9Year <- as.factor(paste(dat$region9, dat$Year, sep = "_"))
  dat$Reg4Year <- as.factor(paste(dat$region4, dat$Year, sep = "_"))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$region9 <- as.factor(dat$region9)
  
  dat <- dat %>%
    group_by(region9) %>%
    mutate(uniqSiteYear = length(unique(SiteYear))) %>%
    filter(uniqSiteYear > 1)
  
  temp <- dat
  
  if(sum(temp$Total) < 20) next
  if(length(unique(temp$SiteID)) < 2) next
  if(length(unique(temp$Year)) < 2) next
  
  mod7 <- try(mod7 <- gam(Total ~ 
                            s(SiteID, bs = "re", k = 5) +
                            s(Reg9Year, bs = "re", k = 5) + 
                            ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                            s(cumdegday, bs = "cr", k = 20),
                          family = nb(theta = NULL, link = "log"),
                          # family = poisson(link = "log"),
                          data = temp,
                          method = "REML", 
                          optimizer = c("outer", "newton"), 
                          sp = c(-1, -1, -1, -1, -1),
                          gamma = 1.4, 
                          control = list(maxit = 500)))
  
  if(("gam" %in% class(mod7)) == FALSE){
    mod7 <- try(mod7 <- gam(Total ~ 
                              s(SiteID, bs = "re", k = 5) +
                              s(Reg9Year, bs = "re", k = 5) + 
                              # ti(Reg9Year, cumdegday, bs = c("re", "cr"), k = c(5, 10)) +
                              s(cumdegday, bs = "cr", k = 20),
                            family = nb(theta = NULL, link = "log"),
                            # family = poisson(link = "log"),
                            data = temp,
                            method = "REML", 
                            optimizer = c("outer", "newton"), 
                            sp = c(-1, -1, -1),
                            gamma = 1.4, 
                            control = list(maxit = 500)))
    
  }
  
  
  if(("gam" %in% class(mod7)) == FALSE) next
  
  gamlist$modgdd <- mod7
  
  
  temp1 <- temp %>% filter(Ordinal > 70 & Ordinal < 320)
  # 
  # start.surv <- sort(unique(temp1$Ordinal))[1]
  # end.surv <- sort(unique(temp1$Ordinal), decreasing = TRUE)[1]
  # allDays <- expand.grid(c(start.surv:end.surv), unique(temp1$SiteYear))
  # names(allDays) <- c("Ordinal", "SiteYear")
  # t <- str_split_fixed(allDays$SiteYear, pattern = "_", 2)
  # allDays$SiteID <- t[,1]
  # allDays$Year <- t[,2]
  
  start.gdd <- sort(unique(temp1$cumdegday))[1]
  end.gdd <- sort(unique(temp1$cumdegday), decreasing = TRUE)[1]
  allgdd <- expand.grid(seq(start.gdd, end.gdd, by = 25), unique(temp1$SiteYear))
  names(allgdd) <- c("cumdegday", "SiteYear")
  t <- str_split_fixed(allgdd$SiteYear, pattern = "_", 2)
  allgdd$SiteID <- t[,1]
  allgdd$Year <- t[,2]
  
  
  newData <- allgdd
  # newData <- allDays
  geo <- temp1 %>% select(SiteID, lat, lon, Reg9Year, region9, Year) %>% distinct()
  newData <- merge(newData, geo, by = c("SiteID", "Year"), all.x = TRUE, all.y = FALSE)
  newData$listlength <- mean(temp1$listlength, na.rm = TRUE)
  newData$temperature <- mean(temp1$temperature, na.rm = TRUE)
  newData$duration <- mean(temp1$duration, na.rm = TRUE)
  
  newData$GAM.pred <- as.vector(predict.gam(mod7, newData, type = "response"))
  # temp$GAM.pred <- predict.gam(mod7, type = "response")
  # newData$GAM.pred <- as.vector(predict.gam(mod7, newData, type = "response",
  #                                 exclude = "ti(Reg9Year,cumdegday"))
  newData <- data.table(newData)
  newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"] 
  
  c <- ggplot(data = newData, aes(x = cumdegday, y = Gamma, group = SiteID, color = SiteID)) +
    geom_line(size = 1, alpha = .5) + 
    theme_bw() + theme(legend.position = "none") +
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
  # print(c)
  d <- ggplot(data = newData, aes(x = cumdegday, y = GAM.pred, group = SiteID, color = SiteID)) + 
    geom_line(size = 1, alpha = .5) +
    theme_bw() + geom_point(data = temp1, aes(x = cumdegday, y = Total))
  d <- d + facet_wrap( ~ Year, scales = "free_y") + theme(legend.position = "none") +
    ggtitle(paste("GAM Predictions", sp, sep = " "))
  # print(d)
  
  
  pdf(paste("CovGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(c)
  dev.off()
  
  pdf(paste("CovGAMPRED", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(d)
  dev.off()
  
  gamlist$mod <- mod7
  gamlist$preds <- newData
  gamlist$counts <- temp1

  saveRDS(gamlist, paste("gam", sp, ".rds", sep = ""))
}


i <- 53
sp <- species[i, "CommonName"]
gamlist <- readRDS(paste("gam", sp, ".rds", sep = ""))





