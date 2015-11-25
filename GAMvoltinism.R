# Plots for targeting species with flexible voltinism
library(lubridate)
library(mgcv)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(gdata)
library(ggplot2)
library(zoo)


setwd("C:/Users/Tyson/Desktop/Box Sync/Ohio/data2012")
data <- read.csv("data.trim.csv", header = TRUE)
data <- data.table(data)
data <- data[, list(SeqID, SiteID.x, SiteDate, Week, Total, CheckListKey, CommonName)]
setnames(data,"SiteID.x","SiteID")
data[, SiteID := as.character(SiteID)]
data[, SiteDate := ymd(as.character(SiteDate))]


region <- data.table(read.csv("site_region.txt", header = TRUE))
region[, SiteID := gsub(" ", "", as.character(SiteID))]
data <- merge(data, region, by = "SiteID")
setkey(data, SeqID)

surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Region", "Week"), with = FALSE])

site_geo <- read.csv("../GIS/OHsites_reconciled.csv", header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := as.character(SiteID)]

ScaleSumTo1 <- function(x){x/sum(x)}

ScaleDennis <- function(x){exp(x)/sum(exp(x))}


ScaledPhenology <- function(counts, yr, r){
  temp <- counts[Year == yr][Region == r]
  # knts <- 10 # try this until it breaks down
    if (dim(temp)[1] < 10){
      knts <- dim(temp)[1]
    }else{
      knts <- 10
    }
  
  if(length(unique(temp$SiteID)) == 1){
    species.model <- gam(Total ~ s(Ordinal, bs = "cr", k = knts), 
                         family = nb(theta = NULL, link = "log"), data = temp, 
                         method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    
    start.surv <- min(temp$Ordinal)
    end.surv <- max(temp$Ordinal)
    
    GAM.pred <- predict.gam(species.model, data.frame(Ordinal = c(start.surv:end.surv)), type = "response")
    output <- data.frame(Year = yr, Region = r, Ordinal = c(start.surv:end.surv), Gamma = ScaleSumTo1(GAM.pred))
  }else{
    
    species.model <- gam(Total ~ SiteID + s(Ordinal, bs = "cr", k = knts), 
                         family = nb(theta = NULL, link = "log"), data = temp, 
                         method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    
    start.surv <- min(temp$Ordinal)
    end.surv <- max(temp$Ordinal)
    maxSite <- temp$SiteID[temp$YearTotal == max(temp$YearTotal)][1] #doesn't really matter which site
    
    GAM.pred <- predict.gam(species.model, data.frame(Ordinal = c(start.surv:end.surv), SiteID = maxSite), type = "response")
    output <- data.frame(Year = yr, Region = r, Ordinal = c(start.surv:end.surv), Gamma = ScaleSumTo1(GAM.pred))
  }
  return(output)
}

SpeciesList <- readRDS("SpeciesList.rds")

region <- data.table(read.csv("site_region.txt", header = TRUE))
region[, SiteID := gsub(" ", "", as.character(SiteID))]


species <- as.character(SpeciesList$CommonName)
species[4] <- "Azures"
species <- species[67:94]

for (sp in species){
  counts <- readRDS(paste("RDSfiles/rawcounts", sp, ".rds", sep = ""))
  counts <- merge(counts, site_geo, by = "SiteID")
  print(sp)
  
  #Filter dataset
  counts <- counts[DaysOut <= 28][Year >= 1998]
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  #   counts[, Ordinal := yday(SiteDate)]
  #   counts[, Year := year(SiteDate)]
  
  #Data restrictions for coming up with regional phenology
  #Site must have >= 3 seen, >= 10 surveys, more than 3 sites in Region for GAM
  datGAM <- counts[YearTotal >= 1]
  datGAM <- datGAM[SurvPerYear >= 5]
  datGAM[, SitesObserved := length(unique(SiteID)), by = list(Year)]
  datGAM <- datGAM[SitesObserved >= 2]
  
  years <- sort(as.numeric(unique(datGAM$Year)))
  regions <- unique(datGAM$Region)
  phen_all <- data.frame()
  for (yr in years){
    for (r in regions){
      if(dim(datGAM[Year == yr][Region == r])[1] == 0) next
      if(dim(datGAM[Year == yr][Region == r][Total > 1])[1] < 3) next
      print(r)
      print(yr)
      phen <- ScaledPhenology(datGAM, yr, r)
      phen_all <- rbind(phen_all, phen)
    }
  }
  
  sums <- datGAM %>%
    group_by(Year, Region) %>%
    summarise(Count = sum(Total))
  
  phen_data <- merge(phen_all, sums, by = c("Year", "Region"))

  c <- ggplot(data = phen_data, aes(x = Ordinal, y = Gamma, group = Region, label = Count)) + 
    geom_line(aes(color = Region), size = .8) 
  c <- c + theme_bw() +
    facet_wrap( ~ Year, ncol = 3, scales = "free_y") + 
    ggtitle(sp) 
  
#   d <- ggplot(data = sums, aes(x = Year, y = Count, group = Region)) + 
#     geom_line(aes(color = Region)) + theme_bw()
#   
  
  pdf(paste("GAMplots/", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(c)
  dev.off()
}
  
  
  
  
  
  
  
  
  
  
  
mod <- gam(Total ~ SiteID + te(Ordinal, lat, k = knts),
           family = nb(theta = NULL, link = "log"), data = temp, 
           method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))

knts <- c(20, 5)
species <- as.character(SpeciesList$CommonName)
species[4] <- "Azures"
# species <- species[70:94]
sp <- species[33]
counts <- readRDS(paste("data2012/RDSfiles/rawcounts", sp, ".rds", sep = ""))

dat <- counts[Year > 1997]

m <- ggplot(dat, aes(x=Ordinal, y = Total, group = Region, color = Region))
# m + geom_density(fill=NA) + facet_grid(. ~ Region)
m + stat_smooth(method = "gam", formula = y ~ s(x, k = 10), se = FALSE,
                family = nb(theta = NULL, link = "log"), 
                gamma = 1.4) + facet_wrap( ~ Year, ncol = 4, scales = "free_y")

