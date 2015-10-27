#sites:rows
#weeks:columns

#separate matrix for each species

#covariate format?
#vector 



library(lubridate)
library(mgcv)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(gdata)
library(ggplot2)
library(zoo)
library(parallel)
library(readr)
library(reshape)
library(reshape2)

##########
#DATA PREP
##########

# setwd("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012")
data <- fread("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/data.trim.csv", header = TRUE)
# data <- data[, list(SeqID, SiteID.x, SiteDate, Week, Total, CheckListKey, CommonName)]
setnames(data,"SiteID.x","SiteID")
data[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]
data[, SiteDate := ymd(as.character(SiteDate))]

dat <- data[year(SiteDate) == 2008]
dat[, `:=` (WeekPerYear = length(unique(Week))), 
          by = list(SiteID)]
dat <- dat[WeekPerYear > 15]

surveys <- distinct(dat[, c("SiteID", "SiteDate", "Week", "SeqID"), with = FALSE])
surveys <- surveys[year(SiteDate) == 2008]


#select top species for dataset
SpeciesNum <- dat %>%
  group_by(CommonName) %>%
  summarise(Present = length(Total)) %>%
  arrange(-Present) %>%
  data.frame()

SpeciesList <- SpeciesNum[-grep("Unidentified", SpeciesNum$CommonName, fixed = TRUE), ]
SpeciesList <- SpeciesList[-which(SpeciesList$CommonName == "None seen this day"), ]
SpeciesList <- filter(SpeciesList, Present > 50)
saveRDS(SpeciesList, "species_expanded.rds")


dat <- dat[CommonName %in% SpeciesList$CommonName]


species <- as.character(SpeciesList$CommonName)
count_array <- array(NA, dim=c(length(unique(surveys$SiteID)), length(unique(surveys$Week)), length(species)))
for (i in 1:length(species)){
  
  sp <- species[i]
  spdat <- dat[CommonName == sp]
  spdat <- unique(spdat)
#   setkey(spdat, SeqID)
  
  #Add zeros to surveys when species not counted during a survey
  
  test <- merge(surveys, spdat, by = c("SiteID", "SiteDate", "Week", "SeqID"), all.x = TRUE)
  counts <- test[, `:=` (CheckListKey = NULL,
                         CommonName = NULL)]
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

  count_array[,,i] <- count_matrix
}

saveRDS(count_array, "count_array_expanded.rds")

#covariates

#some already calculated in OHdetprob.RMD
#what to do about NA's in covariates?
oldcovs <- fread("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/survey.covariates.csv")
covs <- merge(surveys, oldcovs, by = "SeqID", all.x = TRUE)

#Celsius-Fahrenheit issues
covs[mean.temp < 45]$mean.temp <- covs[mean.temp < 45]$mean.temp * 1.8 + 32

#overachieving volunteers going out more than once a week!
#choose first one, not averaging, to match with counts
covs <- covs %>% 
  group_by(SiteID, Week) %>%
  arrange(SiteDate) %>%
  mutate(Duplicate = 1:length(SeqID)) %>%
  filter(Duplicate == 1)


surv <- surveys %>%
  group_by(SiteID, Week) %>%
  arrange(SiteDate) %>%
  summarise(SiteDate = SiteDate[1])

covs <- merge(surv, covs, by = c("SiteID", "Week", "SiteDate"), all.x = TRUE)

# some NA's, not more than 30 for covariates
# just assign them as mean (even though not perfect for seasonal variables)
covs <- data.frame(covs)
for(i in 5:ncol(covs)){
  covs[is.na(covs[,i]), i] <- mean(covs[,i], na.rm = TRUE)
}

covs$Ztemp <- poly(covs$mean.temp, 2)[, 1]
covs$Ztemp2 <- poly(covs$mean.temp, 2)[, 2]
covs$Zwind <- scale(covs$mean.wind)
covs$Zcloud <- scale(covs$mean.cloud)
covs$Zduration <- scale(covs$duration)
covs$Zhour <- scale(covs$start.hour)
covs$Zspecies <- scale(covs$num.species)

# new idea for phi, include ordinal day for time/age standard instead of week
covs$Zjulian <- scale(yday(covs$SiteDate))

#cast covs as matrix, so NA's inserted for missing surveys
cov_array <- array(, dim=c(length(unique(surv$SiteID)), length(unique(surv$Week)), 8))

cov_molten <- melt(covs, id = c("SiteID", "Week", "SiteDate"))
cov_array[,,1] <- as.matrix(cast(cov_molten[cov_molten$variable == "Ztemp", ], SiteID ~ Week, value = "value"))
cov_array[,,2] <- as.matrix(cast(cov_molten[cov_molten$variable == "Ztemp2", ], SiteID ~ Week, value = "value"))
cov_array[,,3] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zwind", ], SiteID ~ Week, value = "value"))
cov_array[,,4] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zcloud", ], SiteID ~ Week, value = "value"))
cov_array[,,5] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zduration", ], SiteID ~ Week, value = "value"))
cov_array[,,6] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zhour", ], SiteID ~ Week, value = "value"))
cov_array[,,7] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zspecies", ], SiteID ~ Week, value = "value"))
cov_array[,,8] <- as.matrix(cast(cov_molten[cov_molten$variable == "Zjulian", ], SiteID ~ Week, value = "value"))


saveRDS(cov_array, "covariates_array_expanded.rds")




#GDD from Dan (or maybe Leslie/Rick Reeves?)
gdd <- fread("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/GddResultsAllSites_1996_2012.csv")
names(gdd)[1] <- "SiteID"
gdd$SiteID <- formatC(gdd$SiteID, width = 3, format = "d", flag = "0")
gdd$Date <- mdy(gdd$full_date)
gdd$Year <- year(gdd$Date)
gdd_summary <- gdd %>%
  filter(todayGDD >= 0) %>%
  group_by(SiteID, Year) %>%
  summarise(YearGDD = max(GDD),
            SpringGDD = max(GDD[ordinalEndDayOfYear < 100], na.rm = TRUE), #negative infinity popping up????
            SummerGDD = max(GDD[ordinalEndDayOfYear < 200], na.rm = TRUE),
            FallGDD = max(GDD[ordinalEndDayOfYear < 300], na.rm = TRUE))

gdd_summary[SpringGDD == -Inf]$SpringGDD <- NA

#spring and yearly GDD have lowest correlation, but still .62.
gdd_covs <- gdd_summary[Year == 2008]

#TODO
#still have problem of volunteers doing >1 survey per week
#merging covariates with all necessary surveys
#adding latitude
#outputting covariates for modeling

sites <- read_csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled.csv")
names(sites)[1] <- "SiteID"
sites$SiteID <- formatC(sites$SiteID, width = 3, format = "d", flag = "0")

cov_sites <- merge(sites, gdd_covs, by = "SiteID", all.x = TRUE)

cov_sites <- cov_sites[which(cov_sites$SiteID %in% unique(surveys$SiteID)), ]

# gdd missing from Catawba Island 103
# use gdd from closest other site
# Turn this into a function!
rowNA <- which(is.na(cov_sites$SpringGDD))
d <- dist(cbind(cov_sites$lat, cov_sites$lon), upper = TRUE)
dists <- as.matrix(d)[rowNA,]
mindist <- min(dists[dists > 0])
rowReplace <- which(dists == mindist)
cov_sites[rowNA, c("Year","YearGDD", "SpringGDD", "SummerGDD", "FallGDD")] <-  cov_sites[rowReplace, c("Year", "YearGDD", "SpringGDD", "SummerGDD", "FallGDD")]


saveRDS(cov_sites, "covariates_sites_expanded.rds")

# 
# 
# 
# 
# cov_survs <- distinct(data[, c("SiteDate", "Week", "SiteID", "Start.POSIX", "StartTemp", 
#                             "EndTemp", "StartWindMPH", "EndWindMPH", "duration"), with = FALSE])
# cov_survs <- cov_survs[year(Start.POSIX) == 2008]
# covs <- cov_survs %>%
#   group_by(Week, SiteID) %>%
#   summarise(StartHour = hour(Start.POSIX[1]),
#             Temperature = StartTemp[1],
#             Wind = StartWindMPH[1],
#             MinutesSurv = duration[1])
# 
# surv <- surveys %>%
#   group_by(Week, SiteID) %>%
#   summarise(SiteDate = SiteDate[1])
# 
# covs <- merge(surv, covs, by = c("SiteID", "Week"), all.x = TRUE)
# 
# covs <- merge(covs, sites[,c(1,3,4)], by = "SiteID")
# 
# #cast covs as matrix, so NA's inserted for missing surveys
# cov_array <- array(, dim=c(length(unique(surveys$SiteID)), length(unique(surveys$Week)), 5))
# 
# cov_molten <- melt(covs, id = c("SiteID", "Week", "SiteDate"))
# cov_array[,,1] <- as.matrix(cast(cov_molten[cov_molten$variable == "StartHour", ], SiteID ~ Week, value = "value"))
# cov_array[,,2] <- as.matrix(cast(cov_molten[cov_molten$variable == "Temperature", ], SiteID ~ Week, value = "value"))
# cov_array[,,3] <- as.matrix(cast(cov_molten[cov_molten$variable == "Wind", ], SiteID ~ Week, value = "value"))
# cov_array[,,4] <- as.matrix(cast(cov_molten[cov_molten$variable == "MinutesSurv", ], SiteID ~ Week, value = "value"))
# cov_array[,,5] <- as.matrix(cast(cov_molten[cov_molten$variable == "lat", ], SiteID ~ Week, value = "value"))
# 
# saveRDS(cov_array, "covariates_array.rds")
# 
# 
# 
# count_array <- readRDS('C:/Users/Tyson/Dropbox/SESYNC ZGroup Summer2015/Ohio/count_array.rds')
# cov_array <- readRDS('C:/Users/Tyson/Dropbox/SESYNC ZGroup Summer2015/Ohio/covariates_array.rds')
# 




