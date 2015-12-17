# use regional GAM to visualize phenology

# try doing it for all years at once

ScaledPhenologyAll<- function(counts){
  temp <- counts
  temp$Year <- as.factor(as.character(temp$Year))
  geo <- unique(temp[, c("SiteID", "lat", "lon"), with = FALSE])
  #   temp[, AvgSitePop := mean(log(Total + 1/30)), by = SiteID]
  knts <- c(20, 5) # try this until it breaks down
  
#   if (max(temp$lat) - min(temp$lat) < 1){
#     mod <- gam(Total ~ SiteID + s(Ordinal, k = knts[1], by = Year),
#                family = nb(theta = NULL, link = "log"), data = temp, 
#                method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
#   }else{
#     
    mod <- gam(Total ~ SiteID + te(Ordinal, lat, k = knts, by = Year),
               family = nb(theta = NULL, link = "log"), data = temp, 
               method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    
    # big data version
#     mod <- bam(Total ~ SiteID + te(Ordinal, lat, k = knts),
#                family = poisson(link = "log"), data = temp, 
#                method = "GCV.Cp", gamma = 1.4, control = list(maxit = 500))
    
    #       
    #           mod <- gam(Total ~ SiteID + te(Ordinal, lat, k = knts),
    #                  family = negbin(theta = prelim$family$getTheta(TRUE), link = "log"), data = temp, 
    #                  method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
  # }
  
  start.surv <- min(temp$Ordinal)
  end.surv <- max(temp$Ordinal)
  
  maxSite <- temp$SiteID[temp$YearTotal == max(temp$YearTotal)][1] #doesn't really matter which site
  
  newData <- expand.grid(unique(temp$SiteID), c(start.surv:end.surv), unique(temp$Year))
  names(newData) <- c("SiteID", "Ordinal", "Year")
  newData <- merge(newData, geo, by = "SiteID")
  
  newData$GAM.pred <- predict.gam(mod, newData, type = "response")
  newData <- data.table(newData)
  newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteID"]
  
  return(newData)
}





species <- as.character(SpeciesList$CommonName)
species[4] <- "Azures"

sp <- species[48]
counts <- readRDS(paste("RDSfiles/rawcounts", sp, ".rds", sep = ""))
counts <- merge(counts, site_geo, by = "SiteID")
print(sp)

#Filter dataset
counts <- counts[DaysOut <= 28]
counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                         YearTotal = sum(Total)), 
                 by = list(SiteID, Year)]
#   counts[, Ordinal := yday(SiteDate)]
#   counts[, Year := year(SiteDate)]

#Data restrictions for coming up with regional phenology
#Site must have >= 3 seen, >= 10 surveys, more than 3 sites in Region for GAM
datGAM <- counts[YearTotal >= 3]
datGAM <- datGAM[SurvPerYear >= 10]
datGAM[, SitesObserved := length(unique(SiteID)), by = list(Year)]
datGAM <- datGAM[SitesObserved >= 5]

# reduce data size to test
# datGAM <- datGAM[Year %in% 2008:2010]
datGAM$Year <- as.factor(as.character(datGAM$Year))

phen <- ScaledPhenologyAll(datGAM)  
dat_mod <- merge(datGAM, phen, by = c("Year", "Ordinal", "SiteID", "lat", "lon"))

# plot GA

c <- ggplot(data = dat_mod, aes(x = Ordinal, y = GAM.pred, group = SiteID)) + geom_line(aes(color = lat), size = .8) 
c <- c + theme_bw() +
  scale_colour_gradient2(name = "Latitude", midpoint = mean(range(dat_mod$lat)), low = "red", mid = "yellow", high = "blue") +
  facet_wrap( ~ Year, ncol = 1, scales = "free_y") + 
  ggtitle(sp)
c
    

d <- ggplot(data = dat_mod, aes(x = Ordinal, y = GAM.pred)) + geom_point(aes(color = lat), size = .8) 
d <- d + theme_bw() +
  scale_colour_gradient2(name = "Latitude", midpoint = mean(range(dat_mod$lat)), low = "red", mid = "yellow", high = "blue") +
  # facet_wrap( ~ Year, ncol = 1, scales = "free_y") + 
  ggtitle(sp)
d

    
    
    
    
  }