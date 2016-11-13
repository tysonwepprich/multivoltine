# diapause model with extra generations
# TWo directions:
# 1. Use GDD and photoperiod from one flight period to predict size/timing of the next later in season
# 2. Use whether partial brood impacts growth rate

library(mgcv)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(readr)
library(geosphere)
library(data.table)

# gdd <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/daymet/growingDD_Daymet.RDS")
# gdd <- gdd[gdd$year >= 1995, ]
# gdd <- as.data.frame(gdd)
# gdd$date <- strptime(paste(gdd$year, gdd$yday, sep = " "), format = "%Y %j")
# gdd$date <- as.POSIXct(gdd$date)

site_geo <- fread("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled.csv", header = TRUE)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(as.numeric(SiteID), width = 3, format = "d", flag = "0")]

# gdd <- gdd %>%
#   mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))
# 
# gdd <- left_join(gdd, site_geo, by = "SiteID")
# gdd$date <- as.Date(gdd$date)
# 
# gdd <- gdd %>% 
#   mutate(photo = daylength(lat = lat, doy = date))
# 
# # remaining gdd left in season?
# gdd <- gdd %>%
#   group_by(SiteID, year) %>%
#   mutate(remaindegday = max(cumdegday) - cumdegday)
# 
# saveRDS(gdd, "data/gdd_covs.rds")

gdd <- readRDS("data/gdd_covs.rds")

sitegdd <- gdd %>%
  group_by(SiteID, year) %>%
  summarise(maxGDD = max(cumdegday)) %>%
  ungroup() %>%
  arrange(maxGDD)

# plots with daylength at different gdd
# compares norhern/southernmost sites for gdd and photoperiod
pltdat <- gdd %>%
  filter(SiteID %in% c("016", "138"))

p <- ggplot(pltdat, aes(x = cumdegday, y = photo, group = year, color = year))
p + geom_line() + facet_wrap(~ SiteID, ncol = 1) + theme_bw()




# 1st try model with diapause decision to split penultimate brood
# no GDD timing between broods, just sum all to estimate last brood size

# for each species
# data needed:
# penultimate brood daily GAM.predictions, dates, photoperiod, sites, coordinates
# penultimate brood abundances, timing (gdd left in season), site, year
# last brood abundance, timing (gdd left in season), site, year
# other winter weather covariates influencing productivity
# subsequent 1st brood abundance


SpeciesList <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012/SpeciesList.rds")
species <- SpeciesList %>% arrange(Present) %>% select(CommonName)

gamlist <- readRDS("gamlist_withdetcov.rds")

splist <- gamlist[[101]]
rm(gamlist)
gc()
sp <- species[101,]
# test with Giant Swallowtail population model
temp <- splist$turning
# temp <- tt
t <- str_split_fixed(temp$SiteYear, pattern = "_", 2)
temp$SiteID <- t[,1]
temp$Year <- t[,2]

temp$date <- as.Date(strptime(paste(temp$Year, temp$Ordinal, sep = " "), format = "%Y %j"))

###plots to see broods


turns <- temp %>%
  group_by(SiteYear) %>%
  summarise(numpeaks = length(which(tptype == "peak")),
            minweight = min(weight, na.rm = TRUE))

par(mfcol = c(3, 2))
hist(turns$numpeaks)
hist(turns$minweight, breaks = 50)
hist(temp$Ordinal[which(temp$tptype == "peak")], breaks = 50)

tempfilt <- temp %>% 
  filter(weight > .05)

turns <- tempfilt %>%
  group_by(SiteYear) %>%
  summarise(numpeaks = length(which(tptype == "peak")),
            minweight = min(weight, na.rm = TRUE))

hist(turns$numpeaks)
hist(turns$minweight, breaks = 50)
hist(tempfilt$Ordinal[which(tempfilt$tptype == "peak")], breaks = 50)

###
temp <- as.tbl(temp)
temp <- temp %>%
  filter(weight > .05)
peaks <- temp %>%
  filter(tptype == "peak") %>%
  group_by(SiteID, Year) %>%
  mutate(broodnum = 1:length(tptype)) %>%
  group_by(SiteID, Year) %>%
  mutate(maxbrood = max(broodnum)) # check and see if any univoltine estimates
temp <- merge(temp, gdd, by = c("SiteID", "date"), all.x = TRUE, all.y = FALSE)
peaks <- merge(peaks, gdd, by = c("SiteID", "date"), all.x = TRUE, all.y = FALSE)
peaks <- peaks %>%
  select(SiteID, date, Ordinal, GAM.pred, SiteYear, aoc, weight, Year,
         broodnum, maxbrood, cumdegday, lat, lon, photo, remaindegday)
#adjust this by species
firstbrood <- 1
typbrood <- 2
extrabrood <- 3
peaks <- peaks %>% filter(maxbrood %in% c(typbrood, extrabrood))
# 1st brood always the subsequent brood for growth rates


results <- list()
for (i in 1:length(unique(peaks$SiteYear))){
  sy <- unique(peaks$SiteYear)[i]
  t <- peaks %>% filter(SiteYear == sy)
  if(nrow(t) < extrabrood){
    newrow <- t[1,]
    newrow$aoc <- 0
    newrow$weight <- 0
    newrow$cumdegday <- NA
    newrow$remaindegday <- NA
    newrow$photo <- NA
    newrow$Ordinal <- NA
    newrow$GAM.pred <- NA
    newrow$broodnum <- nrow(t) + 1
    t <- rbind(t, newrow)
  }
  # t$broodnum <- paste("brood", t$broodnum, sep = "")
  nextbrood <- peaks %>% filter(SiteID == t$SiteID[1],
                               Year == as.character(as.numeric(t$Year[1]) + 1))
  if(nrow(nextbrood) == 0) next
  
  nextdat <- nextbrood %>% filter(broodnum == 1) %>%
    select(SiteID, date, Ordinal, SiteYear, aoc, Year, cumdegday, photo)
  for (j in 1:extrabrood){
    tcol <- t %>% filter(broodnum == j) %>%
      ungroup() %>%
      select(aoc, weight, cumdegday, remaindegday, photo)
    names(tcol) <- paste(names(tcol), j, sep = "_")
    nextdat <- cbind(nextdat, tcol)
  }
  results[[length(results) + 1]] <- nextdat
}
moddat <- data.table::rbindlist(results)
rm(results)
gc()

### plot GAM.preds to check strange sites
pltdat <- splist$preds
pltdat <- newData
pltdat <- pltdat %>% filter(Year == "2001")
pltcount <- splist$counts
pltcount <- pltcount %>% filter(Year == "2001")

p <- ggplot(pltdat, aes(x = Ordinal, y = GAM.pred, group = SiteID))
p + geom_line(aes(color = SiteID)) + 
  geom_point(data = pltcount, aes(x = Ordinal, y = Total, group = SiteID, color = SiteID))







temp <- splist$counts
# temp <- temp %>% filter(temp$Year == "2001")
temp <- temp %>% 
  group_by(SiteID) %>%
  mutate(YearTot = sum(Total),
         YearSurv = length(Total))
temp <- temp %>% filter(YearTot >= 5)

mod7 <- try(gam(Total ~ #s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                  s(SiteYear, bs = "re") + #s(Year, bs = "re") +
                  s(Year, Ordinal, bs = "fs", m = 1) +
                   te(lat, Ordinal, bs = "tp", k = c(5, 15)), 
                family = nb(theta = NULL, link = "log"),
                # family = poisson(link = "log"),
                data = temp,
                method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))


start.surv <- sort(unique(temp$Ordinal))[1]
end.surv <- sort(unique(temp$Ordinal), decreasing = TRUE)[1]
allDays <- expand.grid(c(start.surv:end.surv), unique(temp$SiteYear))
names(allDays) <- c("Ordinal", "SiteYear")
t <- str_split_fixed(allDays$SiteYear, pattern = "_", 2)
allDays$SiteID <- t[,1]
allDays$Year <- t[,2]

newData <- allDays
geo <- temp %>% select(SiteID, lat, lon) %>% distinct()
newData <- merge(newData, geo, by = "SiteID", all.x = TRUE, all.y = FALSE)
newData$listlength <- mean(temp$listlength, na.rm = TRUE)
newData$temperature <- mean(temp$temperature, na.rm = TRUE)
newData$duration <- mean(temp$duration, na.rm = TRUE)

newData$GAM.pred <- predict.gam(mod7, newData, type = "response")
# temp$GAM.pred <- predict.gam(mod7, type = "response")
# newData$GAM.pred <- predict.gam(mod7, newData, type = "response", 
#                                 exclude = c("te(lat,lon,Ordinal"))
newData <- data.table(newData)
newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"] 

temp$GAM.pred <-  predict.gam(mod7, type = "response",
                              exclude = c("s(listlength)",
                                          "s(temperature)",
                                          "s(duration)",
                                          "s(SiteID)"))
temp[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"] 

c <- ggplot(data = newData, aes(x = Ordinal, y = Gamma, group = SiteID, color = SiteID)) +
  geom_line(size = .7) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") 
c
d <- ggplot(data = newData, aes(x = Ordinal, y = GAM.pred, group = SiteID, color = SiteID)) + 
  geom_line(size = 1) +
  theme_bw() #+ geom_point(data = temp, aes(x = Ordinal, y = Total))
d <- d + facet_wrap( ~ Year, scales = "free_y") + theme(legend.position = "none")
d















