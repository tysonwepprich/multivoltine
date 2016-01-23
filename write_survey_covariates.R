setwd("C:/Users/Tyson/REPO/multivoltine")
library(lubridate)
library(plyr)
#read and clean data
data <- read.csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/data.trim.csv", header = TRUE)
names(data)[4] <- "SiteID"
data$SiteDate <- ymd(as.character(data$SiteDate))
data$duration <- as.numeric(difftime(data$End.POSIX, data$Start.POSIX, units = "mins"))

# for every survey use plyr functions to find the following
SpeciesList <- function(CheckListKey, Total){
  sp <- CheckListKey
  if (sp[1] == "A099"){
    species <- 0
    un.id.sp <- 0
    abun <- 0
    un.id.abun <- 0
  } else {
    temp <- grep("A",droplevels(sp))
    un.id.sp <- length(temp)
    un.id.abun <- sum(Total[c(temp)])
    species <- length(sp) - un.id.sp
    abun <- sum(Total) - un.id.abun
  }
  return(c(species, un.id.sp, abun, un.id.abun))
}

#runtime ~ 10min
covariates <- ddply(data, .(SeqID), summarise,
                    mean.temp = mean(c(StartTemp[1], EndTemp[1]), na.rm = TRUE), 
                    mean.cloud = mean(c(StartClouds[1], EndClouds[1]), na.rm = TRUE), 
                    mean.wind = mean(c(StartWindMPH[1], EndWindMPH[1]), na.rm = TRUE),
                    start.hour = hour(Start.POSIX[1]), 
                    ord.date = yday(SiteDate[1]),
                    week = Week[1],
                    duration = duration[1],
                    # recorders = Num.Obs[1],
                    num.species = SpeciesList(CheckListKey, Total)[1],
                    un.id.species = SpeciesList(CheckListKey, Total)[2],
                    abund = SpeciesList(CheckListKey, Total)[3],
                    un.id.abund = SpeciesList(CheckListKey, Total)[4])




#remove implausible values
covariates$mean.temp[which(covariates$mean.temp < 20)] <- NA
covariates$mean.cloud[which(covariates$mean.cloud > 100)] <- NA
covariates$mean.wind[which(covariates$mean.wind > 40)] <- NA
#duration issues: some trips really long, because start time not by 24-hour clock, some start time after end time for negative duration
covariates$duration[which(covariates$duration <= 0)] <- NA
covariates$duration[which(covariates$duration >= 300)] <- NA
#issues with start time just ignored, made NA
covariates$start.hour[which(covariates$start.hour >= 19)] <- NA
covariates$start.hour[which(covariates$start.hour <= 7)] <- NA
# covariates$recorders[which(covariates$recorders > 30)] <- NA


write.csv(covariates, "survey.covariates.csv", row.names = FALSE)
summary(covariates)