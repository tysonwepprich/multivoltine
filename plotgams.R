#plot gam models from newest results
#NOTE: modtime useless, only gives seconds and resets each minute

# path <- "E:/R/ButterflyGAM"
path <- "C:/Users/Tyson/Desktop/GAMs/"
gams <- list.files(path)

# gamsfilter <- grep(pattern = "region.loose", x = gams)
gamsfilter <- grep(pattern = "cov.rev.strict", x = gams)
# gamsfilter <- grep(pattern = "extra.strict", x = gams)

# gams<- gams[grep(pattern = "cov", x = gams)]
gams<- gams[grep(pattern = "strict", x = gams)]


######### basic data needs

source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
library(mclust)
ScaleSumTo1 <- function(x){x/sum(x)}

#######

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
#######


##Notes
#Still missing some models (like mourning cloak extra.strict for example)
#Weird differences in plots using covariate model with covariates set to 0, 
#overemphasizes start and end of season when covariates are low

#try predicting with typical values, or rerun models with covariates standardized to week maybe?
#with extra models, early broods not emphasized enough, might be missed in mclust of regional simulations
#also, might be fine to do the mclust by region4 if counts low in some region9, still account for local adaptation in GDD required
# i like cov models, especially for common species, smaller disk space and computation time
# strict cutoff better to ensure enough data for chap3, but still doesn't account for biased lambdas if surveyed sites have 0 sighted

#look as gam smooth of ordinal on duration/list length/temperature and figure out
#how to fix this, maybe standardized by Site (or region) and week (or cumdegday)

#unknown what cutoff to use for chap2 and chap3 species

#takes a long time to run
gamresults <- list()
for (i in 1:length(gams)){
  gamlist <- readRDS(paste(path, gams[i], sep = "/")) 
  
  gamdf <- gamlist[[1]]
  if (is.na(gamlist[[2]][1])){
    gamdf$AIC <- NA
    gamdf$N <- NA
    gamdf$dev.expl <- NA
    gamdf$modtime <- NULL
  }else{
    gamdf$AIC <- AIC(gamlist[[2]])
    gamdf$N <- summary(gamlist[[2]])$n
    gamdf$dev.expl <- summary(gamlist[[2]])$dev.expl
    gamdf$modtime <- NULL
  }
  gamresults[[i]] <- gamdf
}
gamresdf <- bind_rows(gamresults)
resdf <- gamresdf[is.na(gamresdf$species)==FALSE,]
resdf <- resdf %>% group_by(species) %>% mutate(num = 1:length(model))
resdf$model <- paste(resdf$model, resdf$num, sep = ".")
resdf$CommonName <- NULL
resdf$num <- NULL
resdf$cutoff <- NULL
resdf <- resdf %>% arrange(species, AIC) %>% group_by(species) %>% mutate(delta = AIC - min(AIC))
resdf %>% group_by(model) %>% summarise(meandelta = mean(delta))
saveRDS(resdf, "modelsummaryAIC.rds")

for (i in gamsfilter){
  
  
  #input species model
  gamlist <- readRDS(paste(path, gams[i], sep = "/")) 
  if (is.na(gamlist[[2]])){
    next
  }
  #weird workaround for largest gams without updated region.loose models
  # if("Ordinal" %in% dimnames(attr(gamlist[[2]]$terms, "factors"))[[1]]){
  #   next
  # }
  
  gamdf <- gamlist[[1]]
  gamdf$AIC <- AIC(gamlist[[2]])
  gamdf$N <- summary(gamlist[[2]])$n
  gamdf$dev.expl <- summary(gamlist[[2]])$dev.expl
  gamdf$modtime <- NULL
  
  
  #species data for predictions
  sp <- species <- gamdf$species
  model <- gamdf$model
  cutoff <- gamdf$cutoff
  ###############
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
  dat <- dat %>%
    group_by(region9) %>%
    mutate(uniqSiteYear = length(unique(SiteYear))) %>%
    filter(uniqSiteYear > 1)
  
  temp <- dat
  
  # mod2 <- try(gam(Total ~
  #                  # s(listlength)+
  #                  # s(temperature)+
  #                  # s(duration)+
  #                  s(SiteID, bs = "re", k = 5) +
  #                  s(Reg9Year, cumdegday, bs = "fs", k = 5, m = 1) +
  #                  te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 20), d = c(2, 1)) +
  #                  s(Ordinal, bs = "cr", k = 10),
  #                family = nb(theta = NULL, link = "log"),
  #                # family = poisson(link = "log"),
  #                data = temp,
  #                method = "REML",
  #                optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))

  gammod <- mod.extra
  # gammod <- mod1
  datGAM <- temp
  
  
  #with SiteID
  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday, lat, lon, region9, year, meanGDD) %>%
    filter(SiteYear %in% unique(datGAM$SiteYear)) %>%
    filter(yday >= 50 & yday <= 340) %>%
    dplyr::rename(Ordinal = yday,
                  Year = year)
  pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year")]))
  # pred$listlength <- 0
  # pred$temperature <- 0
  # pred$duration <- 0
  pred$zlistlength <- 0
  pred$ztemperature <- 0
  pred$zduration <- 0
  pred$GAM.pred <- as.vector(predict.gam(gammod, pred, type = "response"))
  
  pred$GAM.pred.basic <- as.vector(predict.gam(gammod, pred, type = "response",
                                               exclude = c("s(zlistlength)",
                                                           "s(ztemperature)",
                                                           "s(zduration)",
                                                           "s(SiteID)")))

  pred <- pred %>%
    group_by(SiteYear) %>%
    mutate(SiteYearGDD = max(cumdegday),
           popindex = sum(GAM.pred),
           Gamma = GAM.pred / sum(GAM.pred))
  pred <- pred %>%
    group_by(region9, Year) %>%
    filter(SiteYear == SiteYear[which(popindex == max(popindex))][1])

  pred <- as.data.frame(pred)
  
  
  
  #region9 gams not that useful for seeing broods
################
  ####################### #region9
  # regs <- unique(datGAM$SiteYear)
  # pred <- gdd %>% 
  #   filter(SiteYear %in% regs) %>% 
  #   group_by(SiteYear, region9) %>% 
  #   summarise(maxgdd = max(cumdegday))
  # pred <- merge(pred, data.frame(cumdegday = seq(0, max(pred$maxgdd), 25)))
  # pred <- pred %>% 
  #   group_by(SiteYear, region9) %>% 
  #   filter(cumdegday <= maxgdd[1]) %>% 
  #   arrange(region9, SiteYear, cumdegday, maxgdd)
  # 
  # pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response",
  #                                               exclude = "s(SiteYear)"))
  # pred <- pred %>%
  #   group_by(region9) %>%
  #   filter(maxgdd == max(maxgdd)) %>% 
  #   mutate(region9GDD = max(cumdegday),
  #          Gamma = GAM.pred.reg9yr / sum(GAM.pred.reg9yr)) %>% 
  #   dplyr::select(cumdegday, region9, GAM.pred.reg9yr, region9GDD, Gamma) %>% 
  #   ungroup() %>% 
  #   dplyr::distinct()
  # pred <- as.data.frame(pred)
  # 
  # #plots by region over all years
  # c <- ggplot(data = pred, aes(x = cumdegday, y = Gamma, group = region9)) +
  #   geom_line(size = 1) + 
  #   theme_bw() + theme(legend.position = "none") +
  #   facet_wrap( ~ region9, scales = "free_y", ncol = 3) + 
  #   ggtitle(paste("Regional Scaled Phenology", sp, sep = " "))
  # # print(c)
  # 
  # d <- ggplot(data = pred, aes(x = cumdegday, y = GAM.pred.reg9yr, group = region9)) +
  #   geom_line(size = 1, alpha = 1) +
  #   theme_bw() + geom_point(data = datGAM, aes(x = cumdegday, y = Total), alpha = .5)
  # d <- d + facet_wrap( ~ region9, scales = "free_y", ncol = 3) + theme(legend.position = "none") +
  #   ggtitle(paste("Regional Counts and GAM predictions", sp, sep = " "))
  # # print(d)
  tempcounts <- datGAM %>% filter(SiteYear %in% unique(pred$SiteYear)) %>% data.frame()
  pred$Year <- as.factor(as.character(pred$Year))
  
  # #plots by Site and Year
  c <- ggplot(data = pred, aes(x = cumdegday, y = GAM.pred, group = region9, color = region9)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + theme(legend.position = "right") +
    # geom_point(data = tempcounts, aes(x = cumdegday, y = Total), alpha = .3)+
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("GAM Predictions", sp, sep = " "))
  # print(c)
  
  cc <- ggplot(data = pred, aes(x = cumdegday, y = GAM.pred.basic, group = region9, color = region9)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + theme(legend.position = "right") +
    geom_point(data = tempcounts, aes(x = cumdegday, y = Total), alpha = .3)+
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("GAM Predictions", sp, sep = " "))
  # print(c)
  
  cb <- ggplot(data = pred, aes(x = cumdegday, y = GAM.pred, group = Year, color = Year)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + theme(legend.position = "right") +
    facet_wrap( ~ region9, scales = "free_y") + ggtitle(paste("GAM Predictions", sp, sep = " "))
  
  d <- ggplot(data = pred, aes(x = cumdegday, y = GAM.pred, group = Year, color = Year)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + theme(legend.position = "right") +
    geom_point(data = tempcounts, aes(x = cumdegday, y = Total), alpha = .3)+
    facet_wrap( ~ region9, scales = "free_y") + ggtitle(paste("GAM Predictions", sp, sep = " "))
  
  #######################
  pdf(paste(sp, "yearCovExcl", ".pdf", sep = ""), width = 13, height = 8)
  print(c)
  dev.off()
  
  pdf(paste(sp, "yearCovExclCounts", ".pdf", sep = ""), width = 13, height = 8)
  print(cc)
  dev.off()

  pdf(paste(sp, "regionCovExcl", ".pdf", sep = ""), width = 10, height = 6)
  print(cb)
  dev.off()
  
  pdf(paste(sp, "regionCovExclCounts", ".pdf", sep = ""), width = 10, height = 6)
  print(d)
  dev.off()
  
}




###############


mod2 <- try(gam(Total ~ 
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

temp$region9 <- as.factor(temp$region9)
mod3 <- try(gam(Total ~ 
                  s(SiteID, bs = "re", k = 5) +
                  # s(region9, bs = "re", k = 5)+
                  # s(Year, bs = "re", k = 5) + 
                  te(region9, Year, cumdegday, bs = c("re", "re", "cr"), k = c(5, 5, 20)) +
                  # s(cumdegday, bs = "cr", k = 20) +
                  s(Ordinal, bs = "cr", k = 10),
                family = nb(theta = NULL, link = "log"),
                # family = poisson(link = "log"),
                data = temp,
                method = "REML", 
                optimizer = c("outer", "newton"), 
                gamma = 1.4, 
                control = list(maxit = 500)))

mod4 <- try(gam(Total ~               
                  s(SiteID, bs = "re", k = 5) +
                  s(Reg9Year, cumdegday, bs = "fs", k = 5, m = 1) +
                  te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 20), d = c(2, 1)) +
                  s(Ordinal, bs = "cr", k = 10),
                family = nb(theta = NULL, link = "log"),
                # family = poisson(link = "log"),
                data = temp,
                method = "REML", 
                optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))

mod5 <- try(gam(Total ~      
                  s(listlength)+s(temperature)+s(duration)+
                  s(SiteID, bs = "re", k = 5) +
                  s(Reg9Year, cumdegday, bs = "fs", k = 5, m = 1) +
                  te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 20), d = c(2, 1)) +
                  s(Ordinal, bs = "cr", k = 10),
                family = nb(theta = NULL, link = "log"),
                # family = poisson(link = "log"),
                data = temp,
                method = "REML", 
                optimizer = c("outer", "newton"), gamma = 1.4, control = list(maxit = 500)))

