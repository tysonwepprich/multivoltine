rm(list=ls())
# list.of.packages <- c("devtools", "parallel", "plyr", "dplyr", "tidyr", 
#                       "readr", "data.table", "mgcv", "lubridate", "lme4",
#                       "MuMIn", "mclust", "stringr", "ggplot2", "geosphere",
#                       "broom", "sjPlot", "gridExtra", "randomForest", 
#                       "forestFloor")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages) > 0) install.packages(new.packages)
ScaleSumTo1 <- function(x){x/sum(x)}



# library(devtools)

library(mclust)
library(plyr)
library(dplyr)
library(mgcv)
library(lubridate)
library(tidyr)
library(stringr)
library(ggplot2)
library(geosphere)
library(data.table)

source('rsquared.glmm.R')

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

site_geo <- read.csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled.csv", header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))


# visualize GAM predictions from species models
# preddat <- pred %>% filter(Year == 2003)
# a <- ggplot(data = preddat, aes(x = cumdegday, y = GAM.pred.reg9yr, group = SiteID, color = SiteID))
# a + geom_line() + facet_wrap(~ Reg9Year)

gdd <- readRDS("C:/Users/Tyson/Desktop/Box Sync/Ohio/daymet/growingDD_Daymet.RDS")
# gdd <- gdd[gdd$year >= 1995, ]
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
#what about expected gdd by photoperiod/site?
gdd_left <- gdd %>% 
  group_by(SiteID, year) %>% 
  mutate(siteyrtotalgdd = max(cumdegday)) %>% 
  ungroup() %>% 
  mutate(actualgddleft = siteyrtotalgdd - cumdegday) %>% 
  group_by(SiteID, yday) %>% 
  mutate(expgddleft = mean(actualgddleft)) %>% 
  ungroup() %>% 
  mutate(gddmismatch = actualgddleft - expgddleft) %>% 
  dplyr::select(SiteID,year,SiteYear, yday, cumdegday, actualgddleft, expgddleft, gddmismatch) %>% 
  dplyr::rename(Ordinal = yday)

test <- gdd_left %>% 
  filter(Ordinal==210) %>% 
  group_by(year) %>% 
  summarise(yrgddmm=mean(gddmismatch),
            yrgddalready=mean(cumdegday))

test2 <- gdd_left %>% 
  filter(Ordinal==210, year>=1995) %>% 
  group_by(SiteID) %>% 
  summarise(sitegddmm=mean(gddmismatch),
            sitegddalready=mean(cumdegday))
  
# rm(gdd)

# i <- 67
fs <- list.files("gamGDDordinal")
#only try a few species first for mixture model brood separation
mvspec <- c(10, 23, 26, 40, 28, 44, 46, 52, 55, 58, 62, 74, 79, 71, 77, 65)
mvbroodmin <- c(2, 1, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 3, 2, 2, 2)
mvbroodmax <- c(3, 2, 3, 4, 2, 3, 2, 2, 2, 3, 3, 3, 4, 3, 3, 3)


# # worried about uncertainty and bad clustering when low populations at sites
# counts <- list()
# for(i in 1:16){
#   fs <- list.files("gamGDDordinal")
#   gamlist <- readRDS(paste("gamGDDordinal/", fs[mvspec[i]], sep = ""))
#   # gammod <- gamlist$mod
#   # pred <- gamlist$preds
#   temp <- gamlist$counts
#   rm(gamlist)
#   temp$filenum <- i
#   counts[[i]] <- temp
# }

# 
# replot phenology
for(i in 1:16){
  fs <- list.files("gamGDDordinal")
  gamlist <- readRDS(paste("gamGDDordinal/", fs[mvspec[i]], sep = ""))
  sp <- gamlist[[1]]
  gammod <- gamlist$modb
  datGAM <- gamlist$counts
  rm(gamlist)

  #
  datGAM <- temp
  gammod <- mod7b

  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday, meanGDD, lat, lon) %>%
    filter(SiteYear %in% unique(datGAM$SiteYear)) %>%
    filter(yday >= 75 & yday <= 320) %>%
    # filter(yday %in% seq(55, 335, 7)) %>%
    dplyr::rename(Ordinal = yday)
  pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year", "Year")]))


  pred$GAM.pred <- as.vector(predict.gam(gammod, pred, type = "response"))

  newData <- data.table(pred)
  newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"]

  # newData <- newData %>%
  #   group_by(Reg9Year) %>%
  #   mutate(RegGDD = mean(unique(meanGDD))) %>%
  #   filter(SiteYear == SiteYear[1])

  c <- ggplot(data = newData, aes(x = Ordinal, y = Gamma, group = SiteID, color = SiteID)) +
    geom_line(size = 1, alpha = .5) +
    # scale_colour_continuous(high="red", low="blue")+
  theme_bw() + theme(legend.position = "none") +
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
  print(c)
  d <- ggplot(data = newData, aes(x = cumdegday, y = Gamma, group = SiteID, color = SiteID)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + theme(legend.position = "none") +
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
  print(d)

  d <- ggplot(data = newData, aes(x = cumdegday, y = GAM.pred, group = SiteID, color = SiteID)) +
    geom_line(size = 1, alpha = .5) +
    theme_bw() + geom_point(data = temp, aes(x = cumdegday, y = Total))
  d <- d + facet_wrap( ~ Year, scales = "free_y") + theme(legend.position = "none") +
    ggtitle(paste("GAM Predictions", sp, sep = " "))
  print(d)

  pdf(paste("OrdinalGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(c)
  dev.off()

  pdf(paste("gddGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(d)
  dev.off()
  }

# # allcounts <- rbindlist(counts)
# # saveRDS(allcounts, "allcounts.rds")
# allcounts <- readRDS("allcounts.rds")
# allcounts$species <- fs[mvspec[allcounts$filenum]]
# allcounts$species <- gsub('gam.gdd.plus.day', '', allcounts$species)
# allcounts$species <- gsub('.rds', '', allcounts$species)
# siteyrcounts <- allcounts %>% 
#   filter(Ordinal > 60 & Ordinal < 330) %>% 
#   group_by(species, SiteYear) %>% 
#   summarise(tots = sum(Total),
#             numsurv = length(Total),
#             numpres = length(which(Total > 0)),
#             startsurv = min(Ordinal),
#             endsurv = max(Ordinal),
#             lengthsurv = max(Ordinal) - min(Ordinal))
# regsites <- allcounts %>% 
#   group_by(Reg4Year) %>% 
#   summarise(numsites = length(unique(SiteID)))


# Estimates brood clusters for each SiteYear, based on simulations from GAM predictions
# changed this to include days outside monitoring to try 
# and ensure broods near boundary counted by mclust
# for (i in 1:length(fs)){
# BroodMixMod <- function(i){
#   fs <- list.files("gamGDDordinal")
#   gamlist <- readRDS(paste("gamGDDordinal/", fs[mvspec[i]], sep = ""))
#   gammod <- gamlist$modb
#   datGAM <- gamlist$counts
#   rm(gamlist)
#   
#   pred <- gdd %>%
#     dplyr::select(SiteID, SiteYear, yday, cumdegday) %>%
#     filter(SiteYear %in% unique(datGAM$SiteYear)) %>% 
#     filter(yday >= 50 & yday <= 340) %>% 
#     dplyr::rename(Ordinal = yday)
#   pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year")]))
#   
#   pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
#   # pred$GAM.pred.reg9yr2 <- as.vector(predict.gam(gammod, pred, type = "response",
#                                                 # exclude = c("s(SiteID)")  ))
#   pred <- as.data.frame(pred)
#   rm(gdd)
#   rm(datGAM)
#   rm(gammod)
#   
#   nsim <- 500
#   outlist <- list()
#   outdf <- data.frame()
#   for(ns in 1:nsim){
#     for (j in 1:length(unique(pred$SiteYear))){
#       index <- unique(pred$SiteYear)[j]
#       clustdat <- pred %>%
#         filter(SiteYear == index) %>%
#         # select(Ordinal, GAM.pred, lat, SiteID) %>%
#         dplyr::select(Ordinal, cumdegday, GAM.pred.reg9yr) %>%
#         rowwise() %>%
#         mutate(GAMpois = rpois(1, GAM.pred.reg9yr[1])) %>%
#         ungroup()
#  
#       if(sum(clustdat$GAMpois) < 10){
#         next
#       }
#       
#       dd <- rep(clustdat$cumdegday, clustdat$GAMpois)
#       daycount <- rep(clustdat$Ordinal, clustdat$GAMpois)
#       dat <- data.frame(dd, daycount)
#       # 
#       # t <-  try(modbs <- mclustBootstrapLRT(dat$dd, model = "E",
#       #                              nboot = 500,
#       #                              verbose = FALSE))
#       mod.dd <- try(Mclust(dat$dd, G=c(mvbroodmin[i]:mvbroodmax[i]), modelNames = "E"), silent = TRUE)
#       mod.day <- try(Mclust(dat$daycount, G=c(mvbroodmin[i]:mvbroodmax[i]), modelNames = "E"), silent = TRUE)
#       if(class(mod.dd) == "try-error") next
#       # dat$class <- as.character(mod$classification)
#       
#       df1 <- data.frame()
#       for (b in 1:length(mod.dd$parameters$mean)){
#         df1 <- rbind(df1, data.frame(model = "degday",
#                                    nsim = ns,
#                                    SiteYear = index,
#                                    brood = b,
#                                    num = length(which(mod.dd$classification == b)),
#                                    weight = mod.dd$parameters$pro[b],
#                                    mu = mod.dd$parameters$mean[b],
#                                    sigma = sqrt(mod.dd$parameters$variance$sigmasq[1]))
#         )
#       }
#       df2 <- data.frame()
#       for (c in 1:length(mod.day$parameters$mean)){
#         df2 <- rbind(df2, data.frame(model = "ordinal",
#                                    nsim = ns,
#                                    SiteYear = index,
#                                    brood = c,
#                                    num = length(which(mod.day$classification == c)),
#                                    weight = mod.day$parameters$pro[c],
#                                    mu = mod.day$parameters$mean[c],
#                                    sigma = sqrt(mod.day$parameters$variance$sigmasq[1]))
#         )
#       }
#       df <- rbind(df1,df2)
#       outlist[[length(outlist) + 1]] <- df
#       # outdf <- rbind(outdf, df)
#     }
#   }
#   outdf <- data.table::rbindlist(outlist)
#   outdf$filenum <- i
#   return(outdf)
# }
# # }
# 
# test <- lapply(X = c(1:2), FUN = BroodMixMod)
# 
# 
# library(parallel)
# # multicore
# system.time({
#   cl <- makeCluster(4)
#   clusterEvalQ(cl, {
#     library(mclust)
#     library(data.table)
#     library(dplyr)
#     library(mgcv)
#   })
#   clusterExport(cl=cl,
#                 varlist=c("mvspec",
#                           "mvbroodmax",
#                           "mvbroodmin",
#                           "gdd"))
# 
#   test <- parLapply(cl, 1:16, BroodMixMod)
#   stopCluster(cl)
# })
# 
# saveRDS(test, file = "BroodMixSigLim3.rds")
# 
# test <- readRDS("BroodMixSigLim.rds")
# # one question:
# # is rpois from GAM prediction equivalent to rbinom from scaled GAM prediction?

####
#new results from sesync cluster with 500 sims
# a <- readRDS("mclust_results.rds")
brood <- readRDS("mclust_broods.rds")
pop <- readRDS("mclust_pops.rds")

# t <- str_split_fixed(brood$Reg9Year, pattern = "_", 2)
# brood$region9 <- t[,1]
# brood$year <- as.factor(t[,2])

a <- brood %>% 
  group_by(species, Reg9Year, nsim) %>% 
  mutate(weight = num/sum(num)) %>% 
  rename(SiteYear = Reg9Year)

## from what i thought was final cluster analysis
## grouped by reg9year, still has many outliers

# a <- readRDS("mclust_results_reg9yr0.rds")
a <- a %>% filter(model == "degday")
# a <- rbindlist(test)
# rm(test)
# a$species <- fs[mvspec[a$filenum]]
# a$species <- gsub('gam', '', a$species)
# a$species <- gsub('.rds', '', a$species)

# flagging poor clustering models
# compound.clusters <- a %>% 
#   group_by(species, SiteYear, nsim) %>%
#   summarise(mudiff = diff(mu))
zero.cluster <- a %>% 
  group_by(species, SiteYear, nsim) %>% 
  filter(min(num) == 0) %>% 
  mutate(mudiff = abs(mu - mu[which(num == 0)]))  
  
# all cases with num = 0 included in compound clusters
# flags <- a %>% 
#   group_by(species, SiteYear, nsim) %>%
#   summarise(mudiff = diff(mu)) %>%
#   # filter(mudiff <= 101.9 |
#   #        mudiff >= 1142.488)
#   filter(mudiff<=20)
a.cut <- anti_join(a, zero.cluster, by = c("species", "SiteYear", "nsim"))
#this only seemed to change Horace's DW, maybe bc
#it has really tight brood 3 and 4
#did not clarify cases where ~1/2 the sims split bt extra or not

# based on the above exploration,
# I will take clusters with 0 estimates, and combine them with nearest cluster
# mergeclusters <- function(SiteYr, spec, simnumb){
#   temp <- a %>% filter(SiteYear == SiteYr[1],
#                        species == spec[1],
#                        nsim == simnumb[1])
#   if(min(temp$num) > 0){
#     return(temp %>% select(-SiteYear, -species, -nsim))
#   }else{
#     temp0 <- temp %>% 
#       filter(num == 0)
#     tempret <- temp %>% 
#       filter(num > 0)
#     ind <- which(abs(tempret$mu - temp0$mu) == min(abs(tempret$mu - temp0$mu)))
#     tempret$weight[ind] <- tempret$weight[ind] + temp0$weight
#     tempret$mu[ind] <- mean(c(tempret$mu[ind], temp0$mu))
#     tempret <- tempret %>% 
#       arrange(mu) %>% 
#       mutate(brood = 1:length(brood))
#     return(tempret %>% select(-SiteYear, -species, -nsim))
#     
#   }
# }

#too slow with mergecluster function
system.time({
  aa <- zero.cluster %>% 
    group_by(species, SiteYear, nsim) %>%
    mutate(mudiffrank = rank(mudiff)) %>% 
    filter(mudiffrank %in% c(1,2))
  bb <- aa %>% 
    group_by(species, SiteYear, nsim) %>%
    mutate(weight_adj = sum(weight),
           mu_adj = (mu[1]*weight[1] + mu[2]*weight[2])/sum(weight))
  cc <- bb %>% 
    filter(num > 0) %>% 
    mutate(weight = weight_adj,
           mu = mu_adj) %>% 
    dplyr::select(-weight_adj, -mu_adj, -mudiffrank)
  dd <- anti_join(zero.cluster, aa)
  ee <- rbind(cc, dd)
  ee <- ee %>% 
    group_by(species, SiteYear, nsim) %>%
    arrange(mu) %>% 
    mutate(brood = 1:length(brood)) %>% 
    arrange(species, SiteYear, nsim, brood) %>% 
    dplyr::select(-mudiff)
})




# system.time({aa <- zero.cluster %>%
#   group_by(species, SiteYear, nsim) %>%
#   do(mergeclusters(.$SiteYear, .$species, .$nsim))
# })
# saveRDS(aa, "Broods_fixed.rds")
# saveRDS(aa, "mclust_results_fixed.rds")
# aa <- readRDS("mclust_results_fixed.rds")

a_edit <- full_join(a.cut, ee)
# a_edit <- left_join(pops_edit, a_edit)
# compound.clusters <- a_edit %>%
#   group_by(species, SiteYear, nsim) %>%
#   mutate(mudiff = c(diff(mu), NA),
#          total = sum(num))
# 
# flags <- compound.clusters %>%
#   filter(mudiff >= 101.9 &
#          mudiff <= 1142.488)
# 
# 
# pops <- compound.clusters %>% 
#   group_by(species, SiteYear) %>%
#   filter(brood == 1) %>% 
#   summarise(meanpop = mean(total))
# 
# # remove lowest 10% population siteyear for each species
# # rough way to try and get sites with better data
# pops_edit <- pops %>%
#   group_by(species) %>%
#   mutate(n = ntile(meanpop, 10)) %>%
#   filter(!n %in% c(1))
# # doesn't do much to affect the following plots
# 
# 
# # still lots of clusters too close (as far as gdd go between broods)
# cplot <- ggplot(compound.clusters, aes(x = mudiff, group = as.factor(brood), color = as.factor(brood))) +
#   geom_density() +
#   facet_wrap(~species, scales = "free_y") +
#   ggtitle("Estimated GDD between brood and next brood for each of 100 simulations per Site x Year")
# cplot
# 
# totplot <- ggplot(pops, aes(x = meanpop)) +
#   geom_density() +
#   facet_wrap(~species, scales = "free")
# totplot
# 
# muplot <- ggplot(a_edit, aes(x = mu, group = as.factor(brood), color = as.factor(brood))) +
#   geom_density() +
#   facet_wrap(~species, scales = "free_y") +
#   ggtitle("Estimated GDD for each brood for each of 100 simulations per site x year")
# muplot
# 
# 

#try using regionwide G
# regions <- site_geo %>%
#   select(SiteID, region9, region4)



# b <- a.cut %>% 
# b <- a_edit %>%
b <- a %>% 
  group_by(species, SiteYear, nsim) %>%
  mutate(numbrood = max(brood))
b <- b %>%
  group_by(species) %>%
  do( complete(., nesting(SiteYear, nsim, numbrood), brood, 
               fill = list(num = 0, weight = 0, mu = NA, sigma = NA))) %>% 
  data.frame()

# t <- str_split_fixed(b$SiteYear, pattern = "_", 2)
# b$SiteID <- t[,1]
# b$Year <- as.factor(t[,2])
# b <- merge(b, regions, by = "SiteID")
# b$reg9yr <- paste(b$region9, b$Year, sep = "_")

c <- b %>%
  group_by(species, SiteYear, nsim) %>%
  summarise(maxbrood = numbrood[1]) %>%
  group_by(species, SiteYear, maxbrood) %>%
  summarise(nselected = n()) %>%
  group_by(species, SiteYear) %>%
  dplyr::mutate(probsel = nselected / sum(nselected)) %>% 
  data.frame()
c <- dplyr::rename(c, numbrood = maxbrood)

# c.reg <- b %>%
#   group_by(species, reg9yr, SiteYear, nsim) %>%
#   summarise(maxbrood = numbrood[1]) %>%
#   group_by(species, reg9yr, maxbrood) %>%
#   summarise(nselected = n()) %>%
#   group_by(species, reg9yr) %>%
#   mutate(probsel = nselected / sum(nselected))
# c.reg <- dplyr::rename(c.reg, numbrood = maxbrood)
# 

d <- merge(b, c, by = c("species", "SiteYear", "numbrood"))
# d <- merge(b, c.reg, by = c("species", "reg9yr", "numbrood"))

e <- d %>%
  group_by(species, SiteYear, brood, numbrood) %>%
  dplyr::summarise(meanN = median(num),
            meanweight = median(weight),
            meanmu = median(mu, na.rm = TRUE),
            meansigma = median(sigma, na.rm = TRUE),
            sdN = sd(num),
            sdweight = sd(weight),
            sdmu = sd(mu, na.rm = TRUE),
            sdsigma = sd(sigma, na.rm = TRUE),
            prob = probsel[1]) %>% 
  data.frame()
e <- e %>% filter(is.na(species) == FALSE)
#sdmu gives NaN when only 1 nsim for a brood

# try to take simulated site pops from mclust-regional
# combine with regional estimates
# pops <- readRDS("mclust_pops.rds")
pops <- readRDS("provisionalN_SiteYear.rds")
t <- str_split_fixed(pops$SiteYear, pattern = "_", 2)
pops$SiteID <- t[,1]
pops$year <- as.factor(t[,2])
pops <- merge(pops, sites[, c("SiteID", "region9")], by = "SiteID")
pops <- pops %>%
  # group_by(species, SiteYear, SiteID,year,region9) %>% 
  # summarise(meanN=mean(SiteN)) %>% 
  # ungroup() %>% 
  mutate(reg9yr = paste(region9, year, sep = "_")) %>% 
  group_by(species, reg9yr) %>%
  mutate(regN = sum(meanN)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(siteNperc = meanN/regN) %>%
  dplyr::select(species, SiteYear, reg9yr, siteNperc)
e1 <- e %>% 
  dplyr::rename(reg9yr = SiteYear)
t <- str_split_fixed(e1$reg9yr, pattern = "_", 2)
e1$region9 <- t[,1]
e1$year <- as.factor(t[,2])
# pops <- pops %>% 
#   dplyr::select(species, SiteYear, reg9yr, meanN) %>% 
#   rename(meansiteN = meanN)
test <- merge(pops, e1, by = c("species", "reg9yr"), all.x = TRUE)
e <- test %>%
  rowwise() %>%
  mutate(meanN = siteNperc*meanN)

# shows mean gdd estimated in "e", but split by numbrood estimated
# ee <- e %>% 
#   group_by(species, SiteYear, numbrood) %>% 
#   mutate(totalN = sum(meanN)) %>% 
#   filter(totalN >= 50)

# weak effects of filtering dataset on prob of choosing numbrood

# ee <- left_join(e, siteyrcounts, by = c("species", "SiteYear"))
# ee <- ee %>% 
#   group_by(species, SiteYear, numbrood) %>% 
#   mutate(totalN = sum(meanN))
#   # filter(tots > 5,
#   #        numsurv > 20,
#   #        numpres > 2,
#   #        startsurv < 120,
#   #        endsurv > 270,
#   #        lengthsurv > 150)
# 
# pp <- ee %>% 
#   group_by(species, SiteYear) %>% 
#   mutate(p = max(prob))
# mod <- glm(p ~ totalN + numsurv + numpres + startsurv +
#              endsurv + lengthsurv, data = pp, family = binomial(link = "logit"))

muplot <- ggplot(e, aes(x = meanmu, group = as.factor(brood), color = as.factor(brood))) +
  geom_density() +
  facet_grid(numbrood~species, scales = "free_y") +
  ggtitle("Estimated GDD for each brood averaged over 100 simulations per site x year")
muplot
# # 
# # above plot includes cases where numbrood prob is low
# # randomly select numbrood based on prob, 1 per siteyear x species
# e.filt <- e %>% 
#   dplyr::filter(brood > 1) %>% 
#   dplyr::mutate(ID = paste(species, SiteYear, numbrood, sep = "_")) %>% 
#   data.frame()
# e.rand <- e %>% 
#   dplyr::group_by(species, SiteYear) %>% 
#   dplyr::filter(brood == 1) %>% 
#   dplyr::sample_n(., size = 1, replace = FALSE, weight = prob) %>% 
#   dplyr::mutate(ID = paste(species, SiteYear, numbrood, sep = "_")) %>% 
#   data.frame()
# 
# e.filt <- e.filt %>% 
#   filter(ID %in% unique(e.rand$ID)) %>% 
#   data.frame()
# e.new <- rbind(e.rand, e.filt)
# e.new <- e.new %>% arrange(species, SiteYear, brood)
# 
# # #regional G selection prob
# # e.reg <- e %>%
# #   group_by(species, SiteYear) %>%
# #   filter(numbrood == max(brood)) %>%
# #   summarise(regprob = prob[1])
# # t <- str_split_fixed(e.reg$SiteYear, pattern = "_", 2)
# # e.reg$SiteID <- t[,1]
# # e.reg$Year <- as.factor(t[,2])
# # e.reg <- merge(e.reg, regions, by = "SiteID")
# # e.reg$Reg9Year <- as.factor(paste(e.reg$region9, e.reg$Year, sep = "_"))
# # e.test <- merge(e.reg, regsites, by = c("Reg9Year"))
# 
# # if e numbrood selected by prob, f needs to be changed
# f <- e.new
# t <- str_split_fixed(f$SiteYear, pattern = "_", 2)
# f$SiteID <- t[,1]
# f$Year <- as.factor(t[,2])
# 
# # for regional cluster
# t <- str_split_fixed(f$SiteYear, pattern = "_", 2)
# f$region9 <- t[,1]
# f$Year <- as.factor(t[,2])
# broods <- merge(f, distinct(sites[, c("region9", "region4")]), by = "region9")
# 
# muplot <- ggplot(e.new, aes(x = meanmu, group = as.factor(brood), color = as.factor(brood))) +
#   geom_density() +
#   facet_wrap(~species, scales = "free_y") +
#   ggtitle("Estimated GDD for each brood averaged over 500 simulations per site x year")
# muplot
# 
# broods <- merge(f, sites, by = "SiteID")
# 
# broods <- broods %>%
#   rename(N = meanN, weight = meanweight) %>% 
#   group_by(species) %>% 
#   mutate(maxbrood = max(brood)) %>% 
#   group_by(species, SiteYear) %>%
#   arrange(brood) %>%
#   mutate(lastprop = N[maxbrood] / (N[maxbrood-1] + N[maxbrood]),
#          lastN = N[maxbrood]) %>% 
#   data.frame()
# 
# 
# #many cases with numbrood<maxbrood give lastprop of INF
# broods <- broods %>% 
#   mutate(lastprop = replace(lastprop, is.finite(lastprop) == FALSE, 0),
#          lastN = replace(lastN, is.finite(lastN) == FALSE, 0)) %>% 
#   data.frame()



# old way, with N, brood weight combined by prob
# using this for new method (regional phenology, site pop size/gdd)
f <- e %>%
  group_by(species, SiteYear, brood) %>%
  summarise(N = sum(meanN * prob),
            weight = sum(meanweight * prob, na.rm = TRUE),
            mu = sum(meanmu * prob)) #this doesn't work for last brood
# keep and use this for penultimate brood though to combine with gddvars


# for siteyr clustering
t <- str_split_fixed(f$SiteYear, pattern = "_", 2)
f$SiteID <- t[,1]
f$Year <- as.factor(t[,2])
broods <- merge(f, sites, by = "SiteID")


# # for regional cluster
# t <- str_split_fixed(f$SiteYear, pattern = "_", 2)
# f$region9 <- t[,1]
# f$Year <- as.factor(t[,2])
# broods <- merge(f, distinct(sites[, c("region9", "region4")]), by = "region9")

broods <- broods %>%
  group_by(species) %>% 
  mutate(maxbrood = max(brood)) %>% 
  group_by(species, SiteYear) %>%
  arrange(brood) %>%
  mutate(lastprop = N[maxbrood] / (N[maxbrood-1] + N[maxbrood]),
         lastN = N[maxbrood]) %>% 
  data.frame()


#few cases with weight = 1 give lastprop of INF
broods <- broods %>% 
  mutate(lastprop = replace(lastprop, is.finite(lastprop) == FALSE, 0),
         lastN = replace(lastN, is.finite(lastN) == FALSE, 0))

saveRDS(broods, file="brooddata.rds")

# plot proportion of last brood by species x region
# specplot <- ggplot(broods, aes(x = lastprop, group = region4, color = region4)) + geom_density() +
#   facet_wrap(~species, scales = "free_y") + theme_bw() +
#   xlab("Proportion of last brood size to sum(last 2 broods' size)")
# specplot


# plot proportion of last brood 
# sitebrood <- broods %>%
#   filter(brood == max(brood)) 
# 
# broodplot <- ggplot(sitebrood, aes(x = lon, y = lat, color = lastprop)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~ Year)
# broodplot

# for siteyear cluster
# tt <- e
# # tt <- e.new
# t <- str_split_fixed(tt$SiteYear, pattern = "_", 2)
# tt$SiteID <- t[,1]
# tt$Year <- as.factor(t[,2])
# rm(t)
# tt <- merge(tt, sites, by = "SiteID")

# alternative for siteyear, did this for new reg phen/site pop method
tt <- merge(f, sites, by = "SiteID")

# # for regional cluster
# # tt <- e
# tt <- e.new
# t <- str_split_fixed(tt$SiteYear, pattern = "_", 2)
# tt$region9 <- t[,1]
# tt$Year <- as.factor(t[,2])
# rm(t)
# # tt <- merge(tt, distinct(sites[, c("region9", "region4")]), by = "region9")


# variable: photoperiod at peaks

# for regional cluster, group gdd vars
# gdd_left <- merge(gdd_left, distinct(sites[, c("SiteYear", "region9")]), by = "SiteYear")
# t <- str_split_fixed(gdd_left$SiteYear, pattern = "_", 2)
# gdd_left$SiteID <- t[,1]
# gdd_left$Year <- as.factor(t[,2])
# rm(t)
# gdd_left <- merge(gdd_left, sites, by = "SiteID")
# gdd_reg <- gdd_left %>% 
#   group_by(region9, Year, Ordinal) %>% 
#   summarise(cumdegday = mean(cumdegday),
#             actualgddleft = mean(actualgddleft),
#             expgddleft = mean(expgddleft),
#             gddmismatch = mean(gddmismatch),
#             lat = mean(lat)) %>% 
#   mutate(SiteYear = paste(region9, Year, sep = "_"))
# gdd2 <- gdd %>% 
#   group_by(region9, year, yday) %>% 
#   summarise(daymeantemp = mean( (tmax..deg.c. + tmin..deg.c.) / 2) ) %>% 
#   mutate(SiteYear = paste(region9, year, sep = "_"))
# 
# # get photoperiod at brood peaks, gdd left in year at peak, temperature around peak
# gdd2varsREG <- function(SiteYr, meanmu){
#   # temp <- gdd_left %>% filter(SiteYear == SiteYr)
#   temp <- gdd_reg %>% filter(SiteYear == SiteYr)
#   
#   ord <- temp$Ordinal[which(temp$cumdegday > meanmu)[1]]
#   lat <- temp$lat[1]
#   photo <- daylength(lat, ord)
#   gddleft <- temp$actualgddleft[which(temp$cumdegday > meanmu)[1]]
#   gddexpected <- temp$expgddleft[which(temp$cumdegday > meanmu)[1]]
#   gddmismatch <- temp$gddmismatch[which(temp$cumdegday > meanmu)[1]]
# 
#   #mean temperature near peak
#   # temp2 <- gdd %>%
#   #   filter(SiteYear == SiteYr) %>%
#   #   mutate(daymeantemp = (tmax..deg.c. + tmin..deg.c.) / 2)
#   temp2 <- gdd2 %>%
#     filter(SiteYear == SiteYr)
#   tempwindow <- 10
#   meantemp <- temp2 %>%
#     filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
#   meantemp <- mean(meantemp$daymeantemp)
# 
#   data_frame(ord,
#              lat,
#              photo,
#              gddleft,
#              gddexpected,
#              gddmismatch,
#              meantemp)
# }
# 
# get photoperiod at brood peaks, gdd left in year at peak, temperature around peak
gdd2vars <- function(SiteYr, meanmu, lat){
  temp <- gdd_left %>% filter(SiteYear == SiteYr)

  ord <- temp$Ordinal[which(temp$cumdegday > meanmu)[1]]
  photo <- daylength(lat, ord)
  gddleft <- temp$actualgddleft[which(temp$cumdegday > meanmu)[1]]
  gddexpected <- temp$expgddleft[which(temp$cumdegday > meanmu)[1]]
  gddmismatch <- temp$gddmismatch[which(temp$cumdegday > meanmu)[1]]

  #mean temperature near peak
  temp2 <- gdd %>%
    filter(SiteYear == SiteYr) %>%
    mutate(daymeantemp = (tmax..deg.c. + tmin..deg.c.) / 2)

  tempwindow <- 5
  meantemp <- temp2 %>%
    filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
  meantemp <- mean(meantemp$daymeantemp)

  data_frame(ord,
             lat,
             photo,
             gddleft,
             gddexpected,
             gddmismatch,
             meantemp)
}
# 
# # ran this already
# system.time({uu <- tt %>%
#   group_by(species, SiteYear, brood, numbrood) %>%
#   do(gdd2vars(.$SiteYear, .$meanmu, .$lat))
# })
# saveRDS(uu, "gddvars.rds")
# # #regional
# # system.time({uu <- tt %>%
# #   group_by(species, SiteYear, brood, numbrood) %>%
# #   do(gdd2varsREG(.$SiteYear, .$meanmu))
# # })
# # saveRDS(uu, "gddvars_reg.rds")

# alternative for siteyear
system.time({uu <- tt %>%
  group_by(species, SiteYear, brood) %>%
  do(gdd2vars(.$SiteYear, .$mu, .$lat))
})
saveRDS(uu, "gddvars5.rds")

uu <- readRDS("gddvars5.rds")
# 
# uu <- uu %>% data.frame()
# vv <- merge(tt, uu, by = c("species", "SiteYear", "brood", "numbrood"))

#alt, mu for last brood WRONG, will filter it
vv <- merge(tt, uu, by = c("species", "SiteYear", "brood"))

saveRDS(vv,"allbroods.rds")
#made choice here to take nsim means of max # broods estimated
#alternative, use regional code below to use e.new, rand selected numbrood
#problem with e.new for Site clusters == crazy zero-inflation
broods.cut <- broods %>%
  dplyr::select(species, SiteYear, brood, N, weight, lastprop, lastN)
vv.cut <- vv %>%
  group_by(species) %>%
  mutate(specmax = max(brood)) %>%
  filter(brood == specmax-1)
moddat <- merge(vv.cut, broods.cut, by = 
                  c("species", "SiteYear", "brood"))

saveRDS(moddat, "officialM1data.rds")
moddat <- readRDS("officialM1data.rds")

# good figure of raw data!!!
# clear to see some trends (and outliers from misclassifying broods)
summdat <- moddat %>% 
  group_by(species, Year, region9) %>% 
  summarise(size = log(sum(N.x)),
            mu.gdd=mean(mu),
            latitude=mean(lat.x),
            ord=mean(ord),
            nsites=length(unique(SiteID)),
            photo=mean(photo),
            temp=mean(meantemp),
            gddmm=mean(gddmismatch),
            lastprop=mean(lastprop))


test <- ggplot(data = summdat, aes(x = ord, 
                                   y = lastprop, 
                                   group = species, 
                                   color = latitude))+
  geom_point()+
  facet_wrap(~species)+
  theme_bw()
test

broods <- readRDS("brooddata.rds")

#shows classifying errors
test2 <- ggplot(data = e, aes(x = meanmu, 
                                   y = meanweight, 
                                   group = species, 
                                   color = as.factor(brood),
                              alpha=prob))+
  geom_point()+
  facet_wrap(~species)+
  theme_bw()
test2
# #regional
# broods.cut <- broods %>%
#   dplyr::select(species, SiteYear, brood, N, weight, lastprop, lastN)
# vv.cut <- vv %>%
#   group_by(species) %>%
#   mutate(specmax = max(brood)) %>%
#   filter(brood == specmax-1)
# moddat <- merge(vv.cut, broods.cut, by = 
#                   c("species", "SiteYear", "brood"))


# 
# #adding this late, get meantemp window from last brood too
# vv.temp <- vv %>%
#   group_by(species) %>%
#   mutate(specmax = max(brood)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   filter(numbrood == specmax)
# vv.temp <- vv.temp %>%
#   filter(brood == numbrood) %>% 
#   mutate(meantempextra = meantemp) %>% 
#   dplyr::select(species, SiteYear, meantempextra)
# moddat <- merge(moddat, vv.temp, by = c("species", "SiteYear"))

# # problem again, with e.new adjustment
# broods.cut <- broods %>%
#   dplyr::select(species, SiteYear, brood, N, weight, lastprop, lastN)
# vv.cut <- vv %>%
#   group_by(species) %>%
#   mutate(specmax = max(brood))
# moddat <- merge(vv.cut, broods.cut, by =
#                   c("species", "SiteYear", "brood"),
#                 all.y = TRUE)
# moddat <- moddat %>%
#   rowwise() %>%
#   filter(brood == (specmax - 1)) %>%
#   data.frame()

#pairs plot shows that ordinal, photoperiod, gddleft/expected 
#are highly correlated, photoperiod has "hook" from before solstice
#temperature around window of peak not
GGally::ggpairs(moddat[, c(23:28, 34)])
GGally::ggpairs(specdat[, c(35, 34, 31, 30)])


modspec <- unique(moddat$species)
library(broom)
library(sjPlot)
library(lme4)
# library(relaimpo)
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# all species together
lmerdat <- moddat %>%
  mutate(zmeanmu = scale_this(mu),
         zphoto = scale_this(photo),
         zlat = scale_this(lat.y),
         ztemp = scale_this(meantemp),
         # ztempextra = scale_this(meantempextra),
         zgddexp = scale_this(gddexpected),
         zgddleft = scale_this(gddleft),
         zgddmm = scale_this(gddmismatch),
         zord = scale_this(ord),
         zlastprop = as.numeric(I(lastprop^(1/3))))
lmerdat$species <- as.factor(lmerdat$species)
lmerdat$region9 <- as.factor(lmerdat$region9)

#within group centering by species
lmerdat <- moddat %>% 
  group_by(species) %>% 
  mutate(zmeanmu = scale_this(mu),
         zphoto = scale_this(photo),
         zlat = scale_this(lat.y),
         ztemp = scale_this(meantemp),
         # ztempextra = scale_this(meantempextra),
         zgddexp = scale_this(gddexpected),
         zgddleft = scale_this(gddleft),
         zgddmm = scale_this(gddmismatch),
         zord = scale_this(ord),
         zlastprop = as.numeric(I(lastprop^(1/3))))
lmerdat$species <- as.factor(lmerdat$species)
lmerdat$region9 <- as.factor(lmerdat$region9)


listcoef <- list()
listmod <- list()
listcoef2 <- list()
listmod2 <- list()
listplot <- list()
listimp <- list()
for (i in 1:length(modspec)){
  # # scaling before species split
  # specdat <- lmerdat %>%
  #   filter(species == modspec[i])

  #scaling after species split
  specdat <- moddat %>%
    ungroup() %>%
    filter(species == modspec[i]) %>%
    mutate(zmeanmu = scale_this(mu),
           zphoto = scale_this(photo),
           zlat = scale_this(lat.y),
           ztemp = scale_this(meantemp),
           # ztempextra = scale_this(meantempextra),
           zgddexp = scale_this(gddexpected),
           zgddleft = scale_this(gddleft),
           zgddmm = scale_this(gddmismatch),
           zord = scale_this(ord),
           zlastprop = I(lastprop^(1/3)))

  try(fit <- lm(zlastprop ~
              (zord+zlat+ztemp+zgddmm)^2,
              data = specdat))
  # try(fit2 <- lmer(zlastprop ~
  #             zord*ztemp +
  #               (1 + zord*ztemp|region9),
  #           data = specdat))
  # 
  # 
  # fit <- glmer(cbind(round(lastN), round(N)) ~
  #                (zord+zgddexp+ztemp)^2 +
  #                ((zord+zgddexp+ztemp)^2|region9), 
  #              data = specdat,
  #              family = binomial(link="logit"),
  #              glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000000))
  # )
  
  outcoef <- tidy(fit)
  outcoef$species <- modspec[i]
  outmod <- glance(fit)
  outmod$species <- modspec[i]
  
  listcoef[[i]] <- outcoef
  listmod[[i]] <- outmod

  listplot[[i]] <- sjp.int(fit,
               show.ci = TRUE,
               type = "eff",
               title = modspec[i],
               mdrt.values = "meansd",
               geom.colors = c("blue", "dark green", "red"),
               prnt.plot = FALSE)
# listimp[[i]] <- car::vif(fit)
  # listimp[[i]] <- calc.relimp(fit)

  # # with random effects
  # outcoef2 <- tidy(fit2)
  # outcoef2$species <- modspec[i]
  # outmod2 <- r.squared.merMod(fit2)
  # outmod2$species <- modspec[i]
  # 
  # listcoef2[[i]] <- outcoef2
  # listmod2[[i]] <- outmod2
}
speccoef <- rbindlist(listcoef)
specmod <- rbindlist(listmod)
# speccoef2 <- rbindlist(listcoef2)
# specmod2 <- rbindlist(listmod2)
# saveRDS(listimp, "M1relaimpo.rds")
# 
# listimp <- readRDS("M1relaimpo.rds")
# # list of variable importance via relaimpo pkg
# outdf <- data.frame()
# for (i in 1:length(listimp)){
#   species <- modspec[i]
#   var.y <- listimp[[i]]@var.y
#   r2 <- listimp[[i]]@R2
#   lmg <- listimp[[i]]@lmg
#   varnames <- names(lmg)
#   avgcoef <- listimp[[i]]@ave.coeffs[, 10]
#   outdf <- rbind(outdf, data.frame(species, var.y, r2, lmg, avgcoef, varnames))
# }
# row.names(outdf) <- NULL
# s <- outdf  %>%  
#   group_by(varnames) %>% 
#   summarise(meanimp = mean(lmg / r2),
#             sdimp = sd(lmg/r2),
#             meanavgcoef = mean(avgcoef),
#             sdavgcoef = sd(avgcoef)) %>% 
#   arrange(-meanimp)
# capture.output(s, file = "M1varimp.txt")


cols <- names(speccoef)[2:5]
speccoef[,(cols) := round(.SD,3), .SDcols=cols]
speccoef %>% data.frame()

s <- speccoef %>% data.frame()
capture.output(s, file = "M1coef.final.txt")
s <- specmod
capture.output(s, file = "M1glance.final.txt")


library(gridExtra)
#pdf of model coefficients
pdf("data_output_specscaled.pdf", height=33, width=8.5)
grid.table(speccoef)
dev.off()


#interaction plots from lm
#3 pdfs: 1 for each variable combo
int.grob <- list()
for (i in 1:16){
  int.grob[[i]] <- listplot[[i]]$plot.list[[1]]
}
n <- length(int.grob)
nCol <- floor(sqrt(n))
pdf("M1interaction1new.pdf", height=15, width=20)
do.call("grid.arrange", c(int.grob, ncol=nCol))
dev.off()

# all species together
lmerdat <- moddat %>%
  mutate(zmeanmu = scale_this(mu),
         zphoto = scale_this(photo),
         zlat = scale_this(lat.y),
         ztemp = scale_this(meantemp),
         # ztempextra = scale_this(meantempextra),
         zgddexp = scale_this(gddexpected),
         zgddleft = scale_this(gddleft),
         zgddmm = scale_this(gddmismatch),
         zord = scale_this(ord),
         zlastprop = I(lastprop^(1/3)))
lmerdat$species <- as.factor(lmerdat$species)
lmerdat$region9 <- as.factor(lmerdat$region9)


#is it ok that scaled predictors have different 
#ranges for different species?
dp <- ggplot(lmerdat, aes(x = ztemp*zord, group = species)) +
  geom_density() +
  facet_wrap(~species)
dp

fit0 <- lmer(zlastprop ~
                  (zord+ztemp+zlat)^2 +
               (1 + (zord+ztemp+zlat)^2|species),
                data = lmerdat)

fit0a <- lmer(zlastprop ~
               (zord+zlat+ztemp+zgddmm)^2 +
               (1 + (zord+zlat+ztemp+zgddmm)^2|species),
             data = lmerdat,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))

null <- gls(zlastprop ~
              (zord+zlat+ztemp+zgddmm)^2,
            data = lmerdat, method = "REML")

fit1 <- lme(zlastprop ~
              (zord+zlat+ztemp+zgddmm)^2,
            data = lmerdat,
            random = ~1|species,
            method = "REML")

# did forward selection for random effects, only one interaction not included
fit1 <- lmer(zlastprop ~
               (zord+zlat+ztemp+zgddmm)^2 +
               (1 
                + zord
                + zlat
                + ztemp
                + zgddmm
                + zord:zlat
                + zord:ztemp
                + zord:zgddmm
                + zlat:ztemp
                + zlat:zgddmm
                + ztemp:zgddmm
                |species),             
             data = lmerdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))
AIC(fit1, fit2)

fit2 <- lmer(zlastprop ~
               (zord+zlat+ztemp+zgddmm)^2 +
               (1 
                + zord
                + zlat
                + ztemp
                + zgddmm
                + zord:zlat
                + zord:ztemp
                + zord:zgddmm
                + zlat:ztemp
                + zlat:zgddmm
                # + ztemp:zgddmm
                |species),             
             data = lmerdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))

AIC(fit1, fit2)


#fit final model with REML
#only one not selected
fit2 <- lmer(zlastprop ~
               (zord+zlat+ztemp+zgddmm)^2 +
               (1 
                + zord
                + zlat
                + ztemp
                + zgddmm
                + zord:zlat
                + zord:ztemp
                + zord:zgddmm
                + zlat:ztemp
                + zlat:zgddmm
                + ztemp:zgddmm
                |species),             
             data = lmerdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))




# NEWEST ATTEMPT AT MODELING
# within species, within region centering
# first scale everything, 
# 1st model: between region variation within each separate species
# filter by species, van der Pol model with region scaling/mean parameters
# 2nd model: between species variation across all regions pooled
# van der Pol model with species scaling/mean parameters
modspec <- unique(moddat$species)
library(broom)
library(sjPlot)
library(lme4)
library(relaimpo)
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


lmerdat <- moddat %>%
  mutate(zmeanmu = scale_this(mu),
         zphoto = scale_this(photo),
         zlat = scale_this(lat.y),
         ztemp = scale_this(meantemp),
         # ztempextra = scale_this(meantempextra),
         zgddexp = scale_this(gddexpected),
         zgddleft = scale_this(gddleft),
         zgddmm = scale_this(gddmismatch),
         zord = scale_this(ord),
         zlastprop = I(lastprop^(1/3)))
lmerdat$species <- as.factor(lmerdat$species)
lmerdat$region9 <- as.factor(lmerdat$region9)


listcoef0 <- list()
listmod0 <- list()
listmodfit0 <- list()
listcoef1 <- list()
listmod1 <- list()
listmodfit1 <- list()
listcoef2 <- list()
listmod2 <- list()
listmodfit2 <- list()
listplot <- list()
for (i in 1:length(modspec)){
  # scaling after species split
  specdat <- moddat %>%
    filter(species == modspec[i]) %>%
    mutate(zmeanmu = scale_this(mu),
           zphoto = scale_this(photo),
           zlat = scale_this(lat.y),
           ztemp = scale_this(meantemp),
           # ztempextra = scale_this(meantempextra),
           zgddexp = scale_this(gddexpected),
           zgddleft = scale_this(gddleft),
           zgddmm = scale_this(gddmismatch),
           zord = scale_this(ord),
           zlastprop = I(lastprop^(1/3))) %>% 
    group_by(region9) %>% 
    mutate(
      var.lat = zlat - mean(zlat),
      var.ord = zord - mean(zord),
      var.temp = ztemp - mean(ztemp),
      var.gddmm = zgddmm - mean(zgddmm),
      mean.ord = mean(zord),
      mean.lat = mean(zlat),
      mean.temp= mean(ztemp),
      mean.gddmm = mean(zgddmm))
  
  #van der pol formula examples
  fit0 <- lmer(zlastprop ~
                 zord+ztemp+zgddmm+
                 (1 
                  # + var.ord
                  # + var.temp
                  # + var.gddmm
                  |region9), 
               data = specdat, REML = TRUE,
               control = lmerControl(optCtrl = list(maxfun = 1e8)))
  
  
  fit1 <- lmer(zlastprop ~
                 var.ord + var.temp + var.gddmm+
                 mean.ord + mean.temp + mean.gddmm+
                 (1 
                  # + var.ord
                  # + var.temp
                  # + var.gddmm
                  |region9), 
               data = specdat, REML = TRUE,
               control = lmerControl(optCtrl = list(maxfun = 1e8)))
  
  fit2 <- lmer(zlastprop ~
                 zord+ztemp+zgddmm+
                 mean.ord + mean.temp+mean.gddmm+
                 (1 
                  # + var.ord
                  # + var.temp
                  # + var.gddmm
                  |region9), 
               data = specdat, REML = TRUE,
               control = lmerControl(optCtrl = list(maxfun = 1e8)))
  
  outcoef <- tidy(fit0)
  outcoef$species <- modspec[i]
  outmodfit <- glance(fit0)
  outmodfit$species <- modspec[i]
  outmodfit <- cbind(outmodfit, r.squared.merMod(fit0))
  
  listcoef0[[i]] <- outcoef
  listmod0[[i]] <- fit0
  listmodfit0[[i]] <- outmodfit
  
  outcoef <- tidy(fit1)
  outcoef$species <- modspec[i]
  outmodfit <- glance(fit1)
  outmodfit$species <- modspec[i]
  outmodfit <- cbind(outmodfit, r.squared.merMod(fit1))
  
  listcoef1[[i]] <- outcoef
  listmod1[[i]] <- fit1
  listmodfit1[[i]] <- outmodfit
  
  outcoef <- tidy(fit2)
  outcoef$species <- modspec[i]
  outmodfit <- glance(fit2)
  outmodfit$species <- modspec[i]
  outmodfit <- cbind(outmodfit, r.squared.merMod(fit2))
  
  listcoef2[[i]] <- outcoef
  listmod2[[i]] <- fit2
  listmodfit2[[i]] <- outmodfit
}

modfits <- rbindlist(listmodfit)
modfits2 <- rbindlist(listmodfit)
modcoef <- rbindlist(listcoef)





# all species model
alldat <- lmerdat %>% 
  group_by(species, region9) %>% 
  mutate(
    var.lat = zlat - mean(zlat),
    var.ord = zord - mean(zord),
    var.temp = ztemp - mean(ztemp),
    var.gddmm = zgddmm - mean(zgddmm),
    mean.ord = mean(zord),
    mean.lat = mean(zlat),
    mean.temp= mean(ztemp),
    mean.gddmm = mean(zgddmm))


# this seems like best model
# no region effect, because mean's cover it
# could use this, or individ. species lm's 
try(allfit <- lmer(zlastprop ~
                  var.ord+var.temp+var.gddmm +
                  mean.ord+mean.temp+mean.gddmm +
                   var.ord:mean.ord +
                  var.temp:mean.temp+
                  var.gddmm:mean.gddmm+
                  (1 
                   + var.ord
                   + var.temp
                   + var.gddmm
                   +mean.ord
                   +mean.temp
                   +mean.gddmm
                   +var.ord:mean.ord
                    + var.temp:mean.temp
                     +var.gddmm:mean.gddmm
                   |species), 
                data = alldat, REML = TRUE,
                control = lmerControl(optCtrl = list(maxfun = 1e8)))
)

# 
# #collinear variables in van der pol method
# #OK solution: reduce window of temp to 10 days rather than 20
# mm <- model.matrix(~var.ord+var.temp+var.gddmm, data = specdat)
# kappa(mm)
# det(cov(mm[,-1]))
# 
# mod <- lm(zlastprop ~
#   var.ord+var.temp+var.gddmm +
#   mean.ord+mean.temp+mean.gddmm, data = specdat)




lmerdat <- lmerdat %>% 
  filter(species == "Silver-spotted Skipper") %>%
  group_by(species, region9) %>% 
  mutate(
         vlat = zlat - mean(zlat),
         vord = zord - mean(zord),
         vtemp = ztemp - mean(ztemp),
         vgddmm = zgddmm - mean(zgddmm),
         meanord = mean(zord),
         meanlat = mean(zlat),
         meantemp= mean(ztemp),
         meangddmm = mean(zgddmm))
lmerdat$species <- as.factor(lmerdat$species)
lmerdat$region9 <- as.factor(lmerdat$region9)

fit1 <- lmer(zlastprop ~
                  vord+vtemp+vgddmm +
                  meanord+meantemp+meangddmm +
                  (1 
                   + vord
                   # + zlat
                   + vtemp
                   + vgddmm
                   + meanord
                   +meantemp
                   +meangddmm
                   # + zord:zlat
                   # + zord:ztemp
                   # + zord:zgddmm
                   # + zlat:ztemp
                   # + zlat:zgddmm
                   # + ztemp:zgddmm
                   |species) ,
                data = lmerdat, REML = TRUE,
                control = lmerControl(optCtrl = list(maxfun = 1e8)))

fit2 <- lmer(zlastprop ~
               vord+vtemp+vgddmm +
               meanord+meantemp+meangddmm +
               (1 
                + vord
                # + zlat
                + vtemp
                + vgddmm
                # +meanord
                # +meantemp
                # +meangddmm 
                # + zord:zlat
                # + zord:ztemp
                # + zord:zgddmm
                # + zlat:ztemp
                # + zlat:zgddmm
                # + ztemp:zgddmm
                |region9), 

             data = lmerdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))

fit5 <- lmer(zlastprop ~
               vord+vtemp+vgddmm +
               meanord+meantemp+meangddmm +
               vord:meanord +vtemp:meantemp + vgddmm:meangddmm +
               (1 
                # + vord
                # # + zlat
                # + vtemp
                # + vgddmm
                # +meanord
                # +meantemp
                # +meangddmm 
                # + zord:zlat
                # + zord:ztemp
                # + zord:zgddmm
                # + zlat:ztemp
                # + zlat:zgddmm
                # + ztemp:zgddmm
                |SiteID), 
             
             data = lmerdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))


fit4 <- lm(zlastprop ~
               vord+vtemp+vgddmm +
               meanord+meantemp+meangddmm +
               vord:meanord +vtemp:meantemp + vgddmm:meangddmm,
             
             data = lmerdat)



AIC(fit1, fit2)


specfit2 <- lmer(zlastprop ~
                   (zord+ztemp+zgddmm+zlat)^2 +
                   meanord+meantemp+meangddmm+meanlat +
                   
                   (1 
                    # + zord
                    # + zlat
                    # + ztemp
                    # + zgddmm
                    # + zord:zlat
                    # + zord:ztemp
                    # + zord:zgddmm
                    # + zlat:ztemp
                    # + zlat:zgddmm
                    # + ztemp:zgddmm
                    |region9),             
                 data = lmerdat, REML = TRUE,
                 control = lmerControl(optCtrl = list(maxfun = 1e8)))




AIC(specfit1, specfit2)


#van der pol formula examples
fit0 <- lmer(zlastprop ~
               zord+ztemp+zgddmm+
               (1 
                # + var.ord
                # + var.temp
                # + var.gddmm
                |region9), 
             data = specdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))


fit1 <- lmer(zlastprop ~
              var.ord + var.temp + var.gddmm+
              mean.ord + mean.temp + mean.gddmm+
              (1 
               # + var.ord
               # + var.temp
               # + var.gddmm
               |region9), 
            data = specdat, REML = TRUE,
            control = lmerControl(optCtrl = list(maxfun = 1e8)))

fit1b <- lm(zlastprop ~
               var.ord + var.temp + var.gddmm+
               mean.ord + mean.temp + mean.gddmm
               , 
             data = specdat)



fit2 <- lmer(zlastprop ~
               zord+ztemp+zgddmm+
               mean.ord + mean.temp+mean.gddmm+
               (1 
                # + var.ord
                # + var.temp
                # + var.gddmm
                |region9), 
             data = specdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))



#convergence issues

library(nloptr)
## from https://github.com/lme4/lme4/issues/98:
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e7)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
  for (n in names(defaultControl)) 
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}
g0.bobyqa2 <- update(fit2,control=lmerControl(optimizer="nloptwrap2"))
g0.NM2 <- update(fitall2,control=lmerControl(optimizer=nloptwrap2,
                                             optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))

s <- summary(fit1)
capture.output(s, file = "M1lmer.final.txt")
saveRDS(fit1, "M1lmer.final.rds")

m1fit <- readRDS("M1lmer.final.rds")
nullfit <- lmer(zlastprop ~
                        1 + (1 
                           + zord
                           + zlat
                           + ztemp
                           + zgddmm
                           + zord:zlat
                           + zord:ztemp
                           + zord:zgddmm
                           + zlat:ztemp
                           + zlat:zgddmm
                           + ztemp:zgddmm
                           |species),             
                        data = lmerdat, REML = TRUE,
                        control = lmerControl(optCtrl = list(maxfun = 1e8)))

sjt.lmer(m1fit, 
         pred.labels = c("Date", "Latitude", "Temperature",
                                "GDD mismatch", "Date x Latitude", "Date x Temperature",
                                "Date x GDD mismatch", "Latitude x Temperature",
                                "Latitude x GDD mismatch", "Temperature x GDD mismatch"),
        depvar.labels = "Response: Proportion of additional brood",
        string.pred = "Predictors",
         p.kr = FALSE,
        show.ci = FALSE, show.se = TRUE,
        show.icc = FALSE, show.re.var = FALSE,
        string.se = "std. error", cell.spacing = .05,
        digits.est = 3, digits.se = 3, file = "m1fit.html")

library(texreg)

htmlreg(m1fit, file = "M1table.doc")

####
#summarise species variation
# with pcagg
library(ggfortify)
coefpca <- coef(m1fit)$species
coefpca$`zlat:zgddmm` <- NULL
coefpca$species <- rownames(coefpca)
df <- coefpca[,c(1:10)]
m1pc <- prcomp(df, center = TRUE, scale. = TRUE)

coefpca <- cbind(coefpca, m1pc$x)

ggplot(data=coefpca, aes(x=PC1,y=PC2,label=species))+
  geom_text()



###phenology

minyr <- lmerdat %>% 
  group_by(species, region9) %>% 
  filter(Year %in% as.character(c(2008, 2010))) %>% 
  summarise(meanbrood = mean(zlastprop), 
            minbrood = min(zlastprop), 
            maxbrood = max(zlastprop)) %>% 
  group_by(species) %>% 
  filter(meanbrood == min(meanbrood))

maxyr <- lmerdat %>% 
  group_by(species, region9) %>%
  filter(Year %in% as.character(c(2008, 2010))) %>% 
  summarise(meanbrood = mean(zlastprop), 
            minbrood = min(zlastprop), 
            maxbrood = max(zlastprop)) %>% 
  group_by(species) %>% 
  filter(meanbrood == max(meanbrood))

# Gonna make a sparktable!
# with species, total obs, site x yr, broods, sparklines

library(sparkTable)

sparkdat <- list()
phendat <- list()
fs <- list.files("gamGDDordinal")
# for (i in 1:16){
for (i in c(1,4,10,13)){
  
  
  gamlist <- readRDS(paste("gamGDDordinal/", fs[sort(mvspec)[i]], sep = ""))
  gammod <- gamlist$mod
  datGAM <- gamlist$counts
  rm(gamlist)
  
  Species <- modspec[i]
  # Total <- sum(datGAM$Total)
  # SiteYears <- length(unique(datGAM$SiteYear))
  # Broods <- paste(mvbroodmin[rank(mvspec)][i], mvbroodmax[rank(mvspec)][i], sep = "-")
  
  
  pred <- gdd %>%
    dplyr::select(SiteID, SiteYear, yday, cumdegday, region4, year, meanGDD) %>%
    filter(SiteYear %in% unique(datGAM$SiteYear)) %>% 
    filter(yday >= 50 & yday <= 340) %>% 
    dplyr::rename(Ordinal = yday)
  pred <- full_join(pred, unique(datGAM[, c("SiteYear", "Reg9Year", "region9")]))
  
  pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
  # pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response",
  #   exclude = c("s(SiteID)")  ))
  pred <- pred %>% 
    group_by(SiteYear) %>% 
    mutate(SiteYearGDD = max(cumdegday),
           Gamma = GAM.pred.reg9yr / sum(GAM.pred.reg9yr))
  pred <- as.data.frame(pred)
  
  # c <- ggplot(data = pred, aes(x = cumdegday, y = Gamma, group = Reg9Year, color = SiteID)) +
  #   geom_line(size = 1, alpha = .5) + 
  #   theme_bw() + theme(legend.position = "none") +
  #   facet_wrap( ~ year, scales = "free_y") #+ ggtitle(paste("Scaled Phenology", Species, sep = " "))
  # print(c)
  
  #output pred for warm/cool years and warm/cool regions for each species
  outregion <- pred %>% 
    # filter(year %in% c(2007, 2008)) %>% 
    group_by(region9, year) %>% 
    summarise(regGDD = mean(unique(SiteYearGDD)),
              uniqSites = length(unique(SiteID)),
              totPop = sum(GAM.pred.reg9yr)) %>%
    ungroup() %>% 
    mutate(Q1 = quantile(totPop, probs = .1),
           Q2 = quantile(totPop, probs = .9)) %>% 
    filter(totPop > Q1) %>%
    group_by(year) %>% 
    filter(uniqSites>2) %>% 
    filter(regGDD == max(regGDD)|regGDD == min(regGDD)) %>% 
    mutate(Reg9Year = paste(region9, year, sep = "_")) %>% 
    group_by(year) %>% 
    mutate(numyr = length(unique(region9)),
           gddyr = mean(regGDD)) %>% 
    ungroup() %>% 
    filter(numyr == 2) %>% 
    filter(gddyr == max(gddyr)|gddyr == min(gddyr)) %>% 
    arrange(gddyr, regGDD) %>% 
    mutate(YearWeather = c("Cool", "Cool", "Warm", "Warm"),
           RegionWeather = c("Cool", "Warm", "Cool", "Warm"))
  outpred <- pred %>% 
    filter(Reg9Year %in% outregion$Reg9Year) %>% 
    # filter(Ordinal <= 315, Ordinal >= 75) %>% 
    group_by(Reg9Year) %>% 
    filter(SiteID == unique(SiteID)[1])
  outpred$species <- Species
  outpred <- merge(outpred, outregion)
  
  phendat[[i]] <- outpred
  # # choose sites with lowest and highest GDD to average phenology
  # cutoffs <- quantile(unique(pred$meanGDD), probs = c(.2, .8))
  #   
  # 
  # cool <- pred %>% 
  #   filter(Ordinal <= 300, Ordinal >= 90, year == 2008) %>% 
  #   mutate(cutoff = quantile(unique(SiteYearGDD), probs = .1)) %>% 
  #   filter(SiteYearGDD < cutoff) %>% 
  #   group_by(SiteYear) %>% 
  #   mutate(GAMpred = GAM.pred.reg9yr / sum(GAM.pred.reg9yr)) %>% 
  #   group_by(Ordinal) %>%
  #   mutate(GAMavg = mean(GAMpred))
  # warm <- pred %>% 
  #   filter(Ordinal <= 300, Ordinal >= 90, year == 2008) %>% 
  #   mutate(cutoff = quantile(unique(SiteYearGDD), probs = .9)) %>% 
  #   filter(SiteYearGDD > cutoff) %>% 
  #   group_by(SiteYear) %>% 
  #   mutate(GAMpred = GAM.pred.reg9yr / sum(GAM.pred.reg9yr)) %>% 
  #   group_by(Ordinal) %>%
  #   mutate(GAMavg = mean(GAMpred))
  # # get fitted values from stat_smooth
  # model <- loess(GAMavg ~ cumdegday, data=cool, span = .15)
  # xrange <- range(cool$cumdegday)
  # xseq <- seq(from=xrange[1], to=xrange[2], length=80)
  # pred <- predict(model, newdata = data.frame(cumdegday = xseq), se=TRUE)
  # y = pred$fit
  # ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
  # ymin = y - ci
  # ymax = y + ci
  # loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
  # fitcool <- loess.DF$y
  # 
  # model <- loess(GAMavg ~ cumdegday, data=warm, span = .15)
  # xrange <- range(warm$cumdegday)
  # xseq <- seq(from=xrange[1], to=xrange[2], length=80)
  # pred <- predict(model, newdata = data.frame(cumdegday = xseq), se=TRUE)
  # y = pred$fit
  # ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
  # ymin = y - ci
  # ymax = y + ci
  # loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
  # fitwarm <- loess.DF$y
  # 
  # 
  # 
  # outdf <- data.frame(Species, Total, SiteYears, Broods, time = c(1:80),  fitcool, fitwarm)
  # sparkdat[[i]] <- outdf
}
# sparkdata <- data.table::rbindlist(sparkdat)
phendata <- data.table::rbindlist(phendat)

t <- str_split_fixed(phendata$Reg9Year, pattern = "_", 2)
phendata$Year <- t[,2]
#rescale gamma
phendata <- phendata %>% 
  group_by(species, YearWeather) %>% 
  mutate(Scaled_Phenology = Gamma / max(Gamma),
         Cumulative_Growing_Degree_Days = cumdegday)

  phendata$YearWeather <- plyr::mapvalues(x = phendata$YearWeather, from = c("Cool", "Warm"),
                  to = c("Cool year", "Warm year"))

# #not quite right, scale off, some species not seen in 2007/08
# quartdata <- phendata %>% 
#   filter(species %in% modspec[1:4])
p <- ggplot(data = phendata, aes(x = Cumulative_Growing_Degree_Days,
                                 y = Scaled_Phenology, color = RegionWeather))
p + geom_line(size = 2) +
  scale_color_manual(values=c("#56B4E9", "#E69F00"))+
  scale_y_continuous(breaks= NULL) +
  scale_x_continuous(breaks = seq(0, 2500, 1000))+
  theme_bw() +
  theme(text = element_text(size=20))+
    # legend.position = "none",
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_blank()) +
  facet_grid(YearWeather ~ species)


phen.grob <- list()
for (i in 1:4){
  quartdata <- phendata %>% 
    filter(species %in% modspec[(4*i-3):(4*i)])
  p <- ggplot(data = quartdata, 
              aes(x = Cumulative_Growing_Degree_Days,
                  y = Scaled_Phenology, color = RegionWeather))+
    geom_line() + 
    scale_color_manual(values=c("#56B4E9", "#E69F00"))+
    scale_y_continuous(breaks= NULL) +
    scale_x_continuous(breaks = seq(0, 2500, 1000))+
    theme_bw() + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    facet_grid(YearWeather ~ species)
  phen.grob[[i]] <- p
}

pdf("phenologyvariation.pdf", height=9, width=9)
do.call("grid.arrange", c(phen.grob, ncol=1))
dev.off()

# sparkdata$time <- 1:80
# prepare content
content <- list(
  function(x){min(x)},
  function(x){max(x)},
  function(x){as.character(x[1])},
  newSparkLine(lineWidth = 2, pointWidth = 0), 
  newSparkLine(lineWidth = 2, pointWidth = 0) 
  )
names(content) <- c("Total", "SiteYears", "Broods", "Cool phenology", "Warm phenology")
# set variables
 vars <- c("Total", "SiteYears", "Broods", "fitcool", "fitwarm")

 # create the sparkTable object
 spd <- data.frame(sparkdata)
spd$Total <- as.numeric(as.character(spd$Total))
spd$SiteYears <- as.numeric(as.character(spd$SiteYears))
spd <- spd[, c("Species", "time","Total", "SiteYears", "Broods", "fitcool", "fitwarm")]
 
 rsSparkdata <- reshapeExt(spd, 
                           idvar = c("Species", "Total", "SiteYears", "Broods"),
                           timeValues = "time",
                           varying  = list(6))
stab <- newSparkTable(dataObj = rsSparkdata, tableContent = content, varType = vars)
# export(stab, outputType = "html", filename = "sparktable.html")
export(stab, outputType="tex", graphNames="o3",filename="t2")


# big species phenology plot











#####
#phylo gls
#####
library(caper)
moddatphy <- moddat
best.tree <- read.tree(file = "../Chap1-Bfly-Landuse-Climate/data/bestTree.txt")
moddatphy$species <- gsub(" ", "", moddatphy$species, fixed = TRUE)
phydat <- comparative.data(data = moddatphy, best.tree, "species")

rm.tips <- phy.tips$tip[-which(phy.tips$tip %in% spec$tip)]
tr.edit <- drop.tip(tr, rm.tips)

# using ape & nlme
library(ape); library(nlme)
# assuming your data are contained in named vectors x & y
y<-y[names(x)]
fit<-gls(zlastprop ~
           (zord+ztemp+zlat)^2,lmerdat,correlation=corPagel(1,best.tree, fixed=FALSE),method="ML")
# assuming your data are in data frame X, with column names
# var1, var2, etc.; for example:
fit<-gls(var1~var2+var3,X,correlation=corPagel(1,tree, fixed=FALSE),method="ML")
######
#lmer plotting 
######
library(merTools)
mod <- readRDS("M1lmer.gddmm.fit.new.rds")

test <- REsim(mod)

# is there a different between coefs in lmer vs species lm?
f <- fixef(mod)
f <- f[-10] #extra fixed effect
r <- ranef(mod)$species
for (i in 1:length(f)){
  r[, i] <- r[, i] + f[i]
}
r <- cbind(r, data.frame("zlat:zgddmm" = fixef(mod)[10]))

# from scratch
# curious about issues with collinear predictors in lmer
# due to species clumping of predictors values

# plot with covariate values drawn across grand range

newdat1 <- expand.grid(species = unique(lmerdat$species),
                       # zord = seq(-2, 2, 0.1),
                       zord = 0,
                       zlat = c(-2, 2),
                       ztemp = 0,
                       # zgddmm = 0,
                       zgddmm = seq(-2, 2, 0.1),
                      zlastprop=0)

fm1 <- m1fit
mm <- model.matrix(terms(fm1),newdat1)
newdat1$zlastprop <- predict(fm1,newdat1,re.form=NA)
newdat1$zlastprop2 <- predict(fm1, newdat1,re.form=NULL)
newdat1$lastprop <- newdat1$zlastprop^3
newdat1$lastprop2 <- newdat1$zlastprop2^3

alldat <- distinct(newdat1[, c("zgddmm", "zlat", "lastprop")])

ggplot(newdat1, aes(x=zgddmm, y=lastprop2, group = interaction(species, zlat), color = zlat)) +
  geom_line(alpha=.5) + 
  geom_line(data=alldat, aes(x=zgddmm, y = lastprop, group = zlat, color = zlat), size = 2) +
  theme_bw() +
  facet_wrap(~zlat)


#plot part 2 lmer results: species variation
preddat <- lmerdat %>% 
  group_by(species) %>% 
  summarise_at(vars(zord,zlat,ztemp,zgddmm),
               funs(mean, sd))


preddata1 <- expand.grid(species = unique(newdat$species),
                        zord = seq(-3.4, 2.5, .05),
                        zlat = seq(-2, 1, .05),
                        ztemp = 0,
                        zgddmm = 0,
                        zlastprop=0)

preddata <- merge(preddat,preddata1)
# only use PC data for each species within 1 SD of its mean
preddata <- preddata %>%
  rowwise() %>% 
  mutate(zordmin = zord_mean - zord_sd,
         zordmax = zord_mean + zord_sd,
         zlatmin = zlat_mean - zlat_sd,
         zlatmax = zlat_mean + zlat_sd) %>% 
  filter(zord >= zordmin && zord <= zordmax) %>% 
  filter(zlat >= zlatmin && zlat <= zlatmax) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  filter(zlat == min(zlat) | zlat == max(zlat)) %>% 
  mutate(latfact = if_else(zlat < zlat_mean, "S", "N"))


fm1 <- m1fit
mm <- model.matrix(terms(fm1), preddata)
preddata1$zlastbrood <- predict(fm1, preddata1, re.form=NA)
preddata$zlastbrood2 <- predict(fm1, preddata, re.form=NULL)
preddata$zlastbrood <- preddata$zlastbrood^3
preddata$zlastbrood2 <- preddata$zlastbrood2^3

preddata$latfact <- as.factor(as.character(preddata$latfact))
levels(preddata$latfact) <- c("Northern latitude", "Southern latitude")
alldat <- distinct(preddata1[, c("zlastbrood", "zord", "zlat")])
alldat <- alldat %>%
  filter(zord >= -2.6 && zord <= 1.65) %>% 
  filter(zlat == -1.5|zlat == 0.5)
alldat$latfact <- as.factor(as.character(alldat$zlat))
levels(alldat$latfact) <- c("Northern latitude", "Southern latitude")

# 
# specdat <- preddata %>% 
#   group_by(species, zlastbrood) %>% 
#   summarise(PC5 = max(PCmax))

# library(ggrepel)
plt <- ggplot(data=preddata, aes(x=zord, y=zlastbrood2, 
                                 group = interaction(species, latfact))) +
  geom_line(alpha=.3, size = .8, color = "black") + 
  facet_grid(~latfact) +
  # geom_text_repel(
  #   data = specdat,
  #   aes(label = species),
  #   size = 6) +
  theme(legend.position = "none") +
  labs(x = "Scaled ordinal date",
       y = "Proportional size of extra brood") +
  geom_line(data=alldat, aes(x=zord, y = zlastbrood, group = latfact), color = "black", size = 1.5) +
  theme_bw()

pdf(paste("M1ordvslat.pdf", sep = ""), width = 4, height = 2)
print(plt)
dev.off()



#merTools
PI.time <- system.time(
  PI <- predictInterval(merMod = fm1, newdata = newdat, 
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
)
preds <- cbind(newdat, PI)

ggplot(preds, aes(x=zgddmm, y=zlastprop2, ymin=lwr, ymax=upr, group = zlat, color = zlat)) +
  geom_point() + 
  geom_linerange() +
  theme_bw() +
  facet_wrap(~species)


sjp.int(mod, type = "eff", facet.grid = TRUE, mdrt.values = "minmax",
        p.kr = FALSE, show.ci = TRUE)


# from glmm wiki
library("lme4")
library("ggplot2") # Plotting
data("Orthodont",package="MEMSS")
fm1 <- lmer(
  formula = distance ~ age*Sex + (age|Subject)
  , data = Orthodont
)
newdat <- expand.grid(
  age=c(8,10,12,14)
  , Sex=c("Female","Male")
  , distance = 0
)
mm <- model.matrix(terms(fm1),newdat)
newdat$distance <- predict(fm1,newdat,re.form=NA)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
tvar1 <- pvar1+VarCorr(fm1)$Subject[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
newdat <- data.frame(
  newdat
  , plo = newdat$distance-cmult*sqrt(pvar1)
  , phi = newdat$distance+cmult*sqrt(pvar1)
  , tlo = newdat$distance-cmult*sqrt(tvar1)
  , thi = newdat$distance+cmult*sqrt(tvar1)
)
#plot confidence
g0 <- ggplot(newdat, aes(x=age, y=distance, colour=Sex))+geom_point()
g0 + geom_errorbar(aes(ymin = plo, ymax = phi))+
  labs(title="CI based on fixed-effects uncertainty ONLY")
#plot prediction
g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
  labs(title="CI based on FE uncertainty + RE variance")




# plot with covariate values drawn across species' ranges











library(effects)
mod <- readRDS("M1lmer.gddmm.fit.new.rds")
ef <- effect("zord", mod)

# dotplot of all random effects
re1 <- ranef(mod, condVar=TRUE, whichel = "species")
lattice::dotplot(re1)

# for each random effects alone
theRan <- ranef(mod, condVar=TRUE)
pv <- attr(theRan$species, "postVar")
se <- pv[1, 1, ]
theIntercepts <- theRan$species[, 1, drop=F]
theFrame <- cbind(theIntercepts, se)
names(theFrame)[1] <- "Intercept"
theFrame$Low <- with(theFrame, Intercept - 2 * se)
theFrame$High <- with(theFrame, Intercept + 2 * se)
theFrame$Variable <- rownames(theFrame)
ggplot(theFrame, 
       aes(y=Intercept, x=reorder(Variable, Intercept))) + 
  geom_linerange(aes(ymin=Low, ymax=High), colour="black") + 
  geom_point(, colour="blue") + 
  coord_flip() + 
  labs(y="Intercept", x=NULL)

#these give neat visualizations of nonlinear effects of
#latitude, temperature, ordinal date
#possible good supplement to linear models.

library(randomForest)
library(forestFloor)

rflist <- list()
for (i in 1:length(modspec)){
  spec <- modspec[i]
  print(spec)
  rfdat <- lmerdat[which(lmerdat$species == spec), ]
  rfdat$zlastprop <- as.numeric(rfdat$zlastprop)
  rfdat <- rfdat %>% 
    dplyr::select(zlastprop, zlat, zord, ztemp, zgddmm)
    # select(zlastprop, meanN, meansigma, meanGDD,
    #        zmeanmu, zphoto, zlat, ztemp, zgddexp, zord)
  y <- rfdat$zlastprop
  x <- rfdat[, -1]
  mod <- randomForest(x = x, y = y, 
                      importance = TRUE, 
                      ntree = 1001, 
                      mtry = 2,
                      keep.inbag = TRUE,
                      keep.forest = TRUE)
  rflist[[i]] <- mod
  rflist[[i]]$species <- spec
}

#compute forestFloor object, often only 5-10% time of growing forest
ff = forestFloor(
  rf.fit = mod,       # mandatory
  X = x,              # mandatory
  calc_np = FALSE,    # TRUE or FALSE both works, makes no difference
  binary_reg = FALSE  # takes no effect here when rfo$type="regression"
)

Col = fcol(ff,cols=1, outlier.lim = 2)

#plot partial functions of most important variables first
plot(ff,                       # forestFloor object
     # plot_seq = ,           # optional sequence of features to plot
     orderByImportance=TRUE,
     plot_GOF = TRUE,
     col = Col
)




#plotting graveyard
summary(fit)
sjp.int(fit, type = "cond")  
r.squared.merMod(fit0)
library("blmeco") 
dispersion_glmer(mod)
library(merTools)
randomSims <- REsim(fit0, n.sims = 500)
# and to plot it
plotREsim(REsim(mod, n.sims = 500), labs = TRUE)
ranks <- expectedRank(mod3, groupFctr = "Year", term = "zphoto")



# messing with climwin
library(climwin)


allbroods <- readRDS(file = "allbroods.rds")
weather <- readRDS(file = "C:/Users/Tyson/REPO/Chap1-Bfly-Landuse-Climate/data/siteweathervars.RDS")
broods <- allbroods %>%
  group_by(species) %>% 
  mutate(maxbrood = max(brood)) %>% 
  group_by(species, SiteYear) %>%
  arrange(brood) %>%
  mutate(lastprop = N[maxbrood] / (N[maxbrood-1] + N[maxbrood]),
         lastN = N[maxbrood]) %>% 
  data.frame()


#few cases with weight = 1 give lastprop of INF
broods <- broods %>% 
  mutate(lastprop = replace(lastprop, is.finite(lastprop) == FALSE, 0),
         lastN = replace(lastN, is.finite(lastN) == FALSE, 0))


newdat <- 
  broods %>%
  group_by(species, SiteYear) %>%
  mutate(siteyearN = sum(N)) %>%
  ungroup %>%
  filter(brood == 1) %>%
  arrange(species, SiteID, Year) %>%
  group_by(species, SiteID) %>%
  mutate(lag.brood1 = lag(N, 1),
         lag.totN = lag(siteyearN, 1),
         lag.siteyear = lag(SiteYear, 1),
         lag.lastprop = lag(lastprop, 1)) %>%
  ungroup() %>%
  mutate(lambda = log((N + 1) / (lag.brood1 + 1))) 

# new error: lambda calculated from nonconsecutive years!!
newdat <- newdat[which(is.na(newdat$lambda)==FALSE),]
lagyr <- str_split_fixed(newdat$lag.siteyear, pattern = "_", 2)
newdat$lagyear <- as.numeric(as.character(lagyr[,2]))
newdat$currentyear <- as.numeric(as.character(newdat$Year))
newdat <- newdat %>% 
  rowwise() %>% 
  mutate(consec = currentyear - lagyear) %>% 
  filter(consec == 1)


# weird errors here with scale
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# break up by species or model together????

#could scale lastprop by site or region, what makes sense?
newdat <- newdat %>%
  group_by(species, SiteID) %>%
  mutate(sitescaleDD1 = scale_this(lag.brood1),
         sitescaleDDall = scale_this(lag.totN),
         specsitescale.prop = scale_this(lag.lastprop^(1/3)))
newdat <- newdat %>%
  group_by(species) %>%
  mutate(specscaleProp = lag.lastprop^(1/3),
         log.lag.brood1 = log(lag.brood1 + 1),
         # log.lag.broodpen = log(lag.penult.brood + 1),
         log.brood1 = log(N + 1))
newdat <- newdat %>% 
  group_by(species, region9) %>% 
  mutate(zlastprop=lag.lastprop-mean(lag.lastprop))

# for sites with too few years, make DD == average or zero
newdat$sitescaleDD1[which(is.na(newdat$sitescaleDD1))] <- 0
newdat$sitescaleDDall <- NULL
newdat$specsitescale.prop <- NULL

# for climwin, need date of lambda
popdat <- newdat %>% 
  dplyr::select(species, SiteYear, SiteID, mu, Year, lat.x, lon,
         region9, ord, zlastprop, lambda, sitescaleDD1) %>% 
  filter(species == "Black Swallowtail")
test <- strptime(paste(popdat$Year, popdat$ord, sep = " "), format = "%Y %j")
popdat$SiteDate <- as.character(as.Date(test, "%d/%m/%Y"), '%d/%m/%Y')


climdat <- gdd %>% 
  mutate(meandailytemp = (tmax..deg.c. + tmin..deg.c.)/2) %>% 
  dplyr::select(SiteID, year, meandailytemp, SiteDate,
                lat,lon,region9, SiteYear)
climdat$SiteDate <- as.character(as.Date(climdat$SiteDate,"%d/%m/%Y"), '%d/%m/%Y')

testmod1 <- slidingwin(
                      exclude = c(2, 2),
                      xvar = list(Temp = climdat$meandailytemp),
                      cdate = climdat$SiteDate,
                      bdate = popdat$SiteDate,
                      baseline = lm(lambda ~ sitescaleDD1, data = popdat),
                      cinterval = "month",
                      range = c(12, 1),
                      type = "absolute",
                      refday = c(1, 6),
                      stat = "mean",
                      func = "quad",
                      # centre = list(popdat$region9, "dev"),
                      cmissing = TRUE,
                      spatial = list(popdat$SiteID, climdat$SiteID)
)

plotdelta(dataset=testmod1[[1]]$Dataset)
plotweights(dataset=testmod1[[1]]$Dataset)
plotbetas(dataset=testmod1[[1]]$Dataset, arrow = TRUE)
plotwin(dataset=testmod1[[1]]$Dataset)

popdat2 <- anti_join(popdat, testmod[[1]]$BestModelData,
                     by = c("lambda"="yvar", 
                            "sitescaleDD1"="sitescaleDD1"))
mod <- lm(lambda~sitescaleDD1 +
          lag.lastprop * (wgdev+wgmean), data =popdat2)

# FLM Teller GAM approach
# aggregating by week gives weird 1/2 week at end of year

clim <- gdd %>% 
  mutate(month = month(SiteDate),
         week = week(SiteDate),
         meanT = (tmax..deg.c.+tmin..deg.c.)/2) %>% 
  dplyr::select(SiteID, year, yday, week, month, meanT, SiteDate, region9, SiteYear)

clim.agg <- clim %>% 
  group_by(SiteID, year, week, region9) %>% 
  summarise(aggT = mean(meanT),
            weekDate = min(SiteDate)) %>% 
  ungroup() %>% 
  group_by(region9, week) %>% 
  mutate(varT = aggT - mean(aggT),
         meanT = mean(aggT))

#lagged climate signals
clim.lag <- clim.agg %>% 
  group_by(SiteID) %>% 
  arrange(SiteID, weekDate) %>% 
  do(data.frame(., setNames(shift(.$aggT, 1:50), 
                            paste("t", formatC(1:50, width=2,flag=0),sep = "."))))
clim.lag <- na.omit(clim.lag)

#when to start climate signals
#based on average species phenology
meanday <- round(mean(popdat$ord))
meanweek <- week(strptime(paste(2015, meanday, sep = " "), format = "%Y %j"))
climdat <- clim.lag %>% 
  filter(week==meanweek) %>% 
  rename(Year=year)

popdat$Year <- as.numeric(as.character(popdat$Year))
alldat <- left_join(popdat, climdat, by = c("SiteID", "Year", "region9"))


## define and enter covariates the way the fitting function wants
tvars <- which(substr(names(alldat),1,2)=="t."); 
alldat$tcovar <- as.matrix(alldat[,tvars]) 

lags <- matrix(0,nrow(alldat),length(tvars)); 
for(i in 1:ncol(lags)) lags[,i]=i; 
alldat$lags=as.matrix(lags); 

lastp <- matrix(0,nrow(alldat),length(tvars)); 
for(i in 1:nrow(lastp)) lastp[i,]=alldat$zlastprop[i]; 
alldat$lastp=as.matrix(lastp); 

knots=40
mod0 <- gam(lambda ~ s(lags, k=knots, by=tcovar, bs="cs"), 
   data=alldat, method="GCV.Cp",gamma=1.4) 
mod1 <- gam(lambda ~ s(sitescaleDD1) + s(lags, k=knots, by=tcovar, bs="cs"), 
            data=alldat, method="GCV.Cp",gamma=1.4) 
mod2 <- gam(lambda ~  s(sitescaleDD1) +
              s(zlastprop) +
              s(lags, k=knots, by=tcovar, bs="cs"), 
            data=alldat, method="GCV.Cp",gamma=1.4) 
mod3 <- gam(lambda ~  s(sitescaleDD1) +
              ti(zlastprop) +
              ti(lags, k=knots, by=tcovar, bs="cs") +
              ti(lastp, lags, by=tcovar), 
            data=alldat, method="GCV.Cp",gamma=1.4) 



a <- ggplot(data = alldat, aes(x=lags, y = tcovar))+geom_line()

#aggregate by month
clim <- gdd %>% 
  mutate(month = month(SiteDate),
         week = week(SiteDate),
         meanT = scale_this((tmax..deg.c.+tmin..deg.c.)/2)) %>% 
  dplyr::select(SiteID, year, yday, week, month, meanT, SiteDate, region9, SiteYear)

clim.agg <- clim %>% 
  group_by(SiteID, year, month, region9) %>% 
  summarise(aggT = mean(meanT),
            monthDate = min(SiteDate)) %>% 
  ungroup() %>% 
  group_by(region9, month) %>% 
  mutate(varT = aggT - mean(aggT),
         meanT = mean(aggT),
         zaggT = as.vector(poly(aggT,2)[,1]),
         zaggT2 = as.vector(poly(aggT,2)[,2]))

#lagged climate signals
clim.lag <- clim.agg %>% 
  group_by(SiteID) %>% 
  arrange(SiteID, monthDate) %>% 
  do(data.frame(., setNames(shift(.$zaggT, 1:12), 
                            paste("t", formatC(1:12, width=2,flag=0),sep = ".")),
                setNames(shift(.$zaggT2, 1:12), 
                         paste("t2", formatC(1:12, width=2,flag=0),sep = "."))))
clim.lag <- na.omit(clim.lag)

#when to start climate signals
#based on average species phenology
meanday <- round(mean(popdat$ord))
meanmonth <- month(strptime(paste(2015, meanday, sep = " "), format = "%Y %j"))
climdat <- clim.lag %>% 
  filter(month==meanmonth) %>% 
  rename(Year=year)

popdat$Year <- as.numeric(as.character(popdat$Year))
alldat <- left_join(popdat, climdat, by = c("SiteID", "Year", "region9"))


## define and enter covariates the way the fitting function wants
tvars <- which(substr(names(alldat),1,2)=="t."); 
alldat$tcovar <- as.matrix(alldat[,tvars]) 

t2vars <- which(substr(names(alldat),1,3)=="t2."); 
alldat$t2covar <- as.matrix(alldat[,t2vars]) 

lags <- matrix(0,nrow(alldat),length(tvars)); 
for(i in 1:ncol(lags)) lags[,i]=i; 
alldat$lags=as.matrix(lags); 

knots=10
mod0 <- gam(lambda ~ s(lags, k=knots, by=tcovar, bs="cs"),
            data=alldat, method="GCV.Cp",gamma=1.4) 
mod1 <- gam(lambda ~ s(lags, k=knots, by=tcovar, bs="cs") +
              s(lags, k=knots, by=t2covar, bs="cs"), 
            data=alldat, method="GCV.Cp",gamma=1.4) 
mod2 <- gam(lambda ~ s(sitescaleDD1) + 
              s(lags, k=knots, by=tcovar, bs="cs") +
              s(lags, k=knots, by=t2covar, bs="cs"), 
            data=alldat, method="GCV.Cp",gamma=1.4) 

alldat$Year <- as.factor(alldat$Year)

# mod2 <- lmer(lambda~
#               sitescaleDD1
#             +zlastprop
#             +t.05
#             # +sitescaleDD1:t.05
#             +zlastprop:t.05
#             +(1 + zlastprop|
#                 Year),
#             data = alldat,REML = TRUE,
#             control = lmerControl(optCtrl = list(maxfun = 1e8)))
#             # control=lmerControl(optimizer="nloptwrap2"))
#               



# original part 2 analysis
# lambda vs climate variables * lastbrood
broods <- readRDS(file = "brooddata.rds")
weather <- readRDS(file = "C:/Users/Tyson/REPO/Chap1-Bfly-Landuse-Climate/data/siteweathervars.RDS")

weather <- data.table(weather)
setnames(weather,"site","SiteID")
setnames(weather, "year", "Year")
weather[, SiteID := formatC(as.numeric(as.character(SiteID)), width = 3, format = "d", flag = "0")]
winter <- weather %>% dplyr::select(SiteID, Year, 
                                    winter_meanTemp, 
                                    # winter_totPrecip, 
                                    fall_meanTemp, 
                                    # fall_totPrecip, 
                                    prevsum_meanTemp,
                                    # prevsum_totPrecip,
                                    # prevspr_meanTemp, 
                                    # prevspr_totPrecip, 
                                    spring_meanTemp)
                                    # spring_totPrecip
                                    # currsum_meanTemp,
                                    # currsum_totPrecip
                                    # )

# use gdd to aggregate temperature by month
# shift dates so easier to summarise over new years
# yearshift corresponds to current year in lambda
test <- gdd %>% 
  mutate(week=week(SiteDate),
         dateshift = SiteDate + duration(6, "months"),
         ydayshift = yday(SiteDate + duration(6, "months")),
         yearshift = year(SiteDate + duration(6, "months")))
# a <- ggplot(test, aes(x = yday, y = numfreeze, group = year, color = year)) + geom_line()
# a
# b <- ggplot(test, aes(x = ydayshift, y = tmin..deg.c., group = year, color = year)) + geom_line()
# b + facet_wrap(~year) + geom_smooth(method="gam", formula = y ~ s(x)) + geom_hline(yintercept = 0)

wintergam <- function(gdd.data){
  g <- gam(tmin..deg.c.~s(ydayshift), data = gdd.data)
  predg <- predict(g, newdata = data.frame(ydayshift = 1:365))
  subzeros <- which(predg < 0)
  winter.onset <- min(subzeros)
  winter.end <- max(subzeros)
  winter.length <- (winter.end - winter.onset)
  winter.min <- min(predg[subzeros])
  winter.cummin <- sum(predg[subzeros])
  temp <- gdd.data %>% filter(ydayshift %in% subzeros)
  winter.snowtot <- sum(temp$swe..kg.m.2.)
  winter.prectot <- sum(temp$prcp..mm.day.)
  return(data.frame(winter.onset, winter.end, winter.length, 
                    winter.min, winter.cummin, winter.snowtot, winter.prectot))
}
winter <- test %>% filter(yearshift >= 1996) %>% 
  group_by(yearshift, SiteID) %>% 
  do(wintergam(.))
#careful with last year, only includes 1/2 winter
winter <- winter %>% 
  filter(yearshift != 2015)

# how correlated are these winter variables?
pc <- prcomp(winter[,-c(1:2)], center = TRUE, scale. = TRUE)
# pc2 <- princomp(winter[,-c(1:2)], cor=TRUE) #this is equivalent
winterdat <- cbind(data.frame(winter), data.frame(pc$x)[,1:5])

# within-site centering on pc's now
# to make sure they include all site x year combos before species added
winterdat <- winterdat %>% 
  group_by(SiteID) %>% 
  mutate(site.var.pc1 = PC1 - mean(PC1),
         site.var.pc2 = PC2 - mean(PC2),
         site.var.pc3 = PC3 - mean(PC3),
         site.var.pc4 = PC4 - mean(PC4),
         site.var.pc5 = PC5 - mean(PC5),
         site.mean.pc1 = mean(PC1),
         site.mean.pc2 = mean(PC2),
         site.mean.pc3 = mean(PC3),
         site.mean.pc4 = mean(PC4),
         site.mean.pc5 = mean(PC5)) %>% 
  ungroup() %>% 
  mutate(site.var.pc1 = scale_this(site.var.pc1),
         site.var.pc2 = scale_this(site.var.pc2),
         site.var.pc3 = scale_this(site.var.pc3),
         site.var.pc4 = scale_this(site.var.pc4),
         site.var.pc5 = scale_this(site.var.pc5),
         site.mean.pc1 = scale_this(site.mean.pc1),
         site.mean.pc2 = scale_this(site.mean.pc2),
         site.mean.pc3 = scale_this(site.mean.pc3),
         site.mean.pc4 = scale_this(site.mean.pc4),
         site.mean.pc5 = scale_this(site.mean.pc5))


# 
# 
# climdat <- gdd %>% 
#   mutate(meandailytemp = (tmax..deg.c. + tmin..deg.c.)/2) %>% 
#   dplyr::select(SiteID, year, meandailytemp, SiteDate,
#                 lat,lon,region9, SiteYear) %>% 
#   mutate(month = month(SiteDate)) %>% 
#   group_by(year, month, SiteID, region9, SiteYear) %>% 
#   summarise(monthmean = mean(meandailytemp))
#   
# climlag <- climdat %>% 
#   group_by(SiteID) %>% 
#   arrange(SiteID, year, month) %>% 
#   do(data.frame(., setNames(shift(.$monthmean, 1:12), 
#                             paste("t", formatC(1:12, width=2,flag=0),sep = "."))))
#               
# climlag <- na.omit(climlag)
# 
# pcdat <- climlag %>% 
#   filter(month == 6,
#          year >= 1996) %>% 
#   data.frame()
# 
# pc <- prcomp(pcdat[,c(7:18)], center = TRUE, scale. = TRUE)
# pcdat <- cbind(pcdat, pc$x)
# 
# 
# a <- ggplot(pcdat, aes(x = year, y = PC3, group = region9))+
#   geom_point()+
#   facet_wrap(~region9)
# a
# model population growth from 1st to 1st broods
# can cite Bradshaw for using this lambda, 
# but also affected by summer conditions that ~last brood prop

# expanddat <- broods %>%
#   tidyr::expand(SiteID, Year) %>%
#   dplyr::left_join(broods)

newdat <- 
  broods %>%
  group_by(species, SiteYear) %>%
  mutate(siteyearN = sum(N)) %>%
  ungroup %>%
  filter(brood == 1) %>%
  arrange(species, SiteID, Year) %>%
  group_by(species, SiteID) %>%
  mutate(lag.brood1 = lag(N, 1),
         lag.totN = lag(siteyearN, 1),
         lag.siteyear = lag(SiteYear, 1),
         lag.lastprop = lag(lastprop, 1)) %>%
  ungroup() %>%
  mutate(lambda = log((N + 1) / (lag.brood1 + 1))) 

# new error: lambda calculated from nonconsecutive years!!
newdat <- newdat[which(is.na(newdat$lambda)==FALSE),]
lagyr <- str_split_fixed(newdat$lag.siteyear, pattern = "_", 2)
newdat$lagyear <- as.numeric(as.character(lagyr[,2]))
newdat$currentyear <- as.numeric(as.character(newdat$Year))
newdat <- newdat %>% 
  rowwise() %>% 
  mutate(consec = currentyear - lagyear) %>% 
  filter(consec == 1)

# try with middle broods
# newdat <- newdat %>% 
#   dplyr::select(species, SiteYear, brood, N, lag.brood1, lag.lastprop)
# newdat2 <- 
#   broods %>%
#   filter(brood == maxbrood-1) %>%
#   arrange(species, SiteID, Year) %>%
#   group_by(species, SiteID) %>%
#   mutate(lag.penult.brood = lag(N, 1),
#          lag.siteyear = lag(SiteYear, 1),
#          lag.lastprop = lag(lastprop, 1),
#          Npen = N) %>% 
#   dplyr::select(species, SiteYear, Npen, lag.penult.brood, lag.siteyear, lag.lastprop, lat, region9, Year)
# newdat <- merge(newdat, newdat2, by = c("species", "SiteYear", "lag.lastprop"))


# weird errors here with scale
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# break up by species or model together????

#could scale lastprop by site or region, what makes sense?
newdat <- newdat %>%
  group_by(species, SiteID) %>%
  mutate(sitescaleDD1 = scale_this(log(lag.brood1+1)),
         # sitescaleDDall = scale_this(lag.totN),
         specsitescale.prop = scale_this(lag.lastprop^(1/3)))
newdat <- newdat %>%
  group_by(species) %>%
  mutate(specscaleProp = lag.lastprop^(1/3),
         log.lag.brood1 = log(lag.brood1 + 1),
         # log.lag.broodpen = log(lag.penult.brood + 1),
         log.brood1 = log(N + 1))

# for sites with too few years, make DD == average or zero
newdat$sitescaleDD1[which(is.na(newdat$sitescaleDD1))] <- 0
newdat$sitescaleDDall <- NULL
newdat$specsitescale.prop <- NULL
# newdat <- newdat %>%
#   group_by(species, region9) %>%
#   mutate(
#          specregscale.prop = scale_this(lag.lastprop^(1/3)))

#set NA to 0 for DD
# newdat$sitescaleDD[is.na(newdat$sitescaleDD)] <- 0
# # newdat$scaleProp[is.na(newdat$scaleProp)] <- 0
# # scaling lastprop
# # from http://stats.stackexchange.com/questions/109702/empirical-logit-transformation-on-percentage-data
# eps <- min(newdat$lastprop[newdat$lastprop > 0], na.rm = TRUE) / 2
# newdat$lastprop[newdat$lastprop == 0] <- eps
# newdat$lastprop[newdat$lastprop == 1] <- 1 - eps
# # doesn't addressed nonnormal, way negative scaled zeros


# newdat$specscalePropSQ <- newdat$specscaleProp ^ 2
# newdat$scalelat <- scale_this(newdat$lat)
# vv.gdd <- vv.cut %>%
#   dplyr::select(species, SiteYear, gddleft, gddmismatch) %>%
#   dplyr::rename(lag.siteyear = SiteYear)
# newdat <- merge(newdat, vv.gdd, by = c("species", "lag.siteyear"))
# newdat$scalegddleft <- scale_this(newdat$gddleft)
# newdat$scalegddmismatch <- scale_this(newdat$gddmismatch)
newdat$region9 <- as.factor(newdat$region9)

# add winter weather 
# winter$Year <- as.character(winter$Year)
# Losing 100 rows of data here????
# should I add prev summer temperature as alternative to brood?
# newdat <- merge(newdat, winter, by = c("SiteID", "Year"))
# newdat$scalewinterT <- scale_this(newdat$winter_meanTemp)
# newdat$scalefallT <- scale_this(newdat$fall_meanTemp)
# newdat$scaleprevsum <- scale_this(newdat$prevsum_meanTemp)
# newdat <- newdat[complete.cases(newdat),]
# 
# # trying PCA with all sites/yrs, THEN merge with newdat
# pcdat <- winter[,-c(1:2), with=FALSE]
# pc <- prcomp(pcdat, center = TRUE, scale. = TRUE)
# pcdat <- cbind(pcdat, pc$x)
# # newpcdat <- merge(winter, pcdat, by = c("winter_meanTemp",
#                                         # "fall_meanTemp",
#                                         # "prevsum_meanTemp",
#                                         # "prevspr_meanTemp",
#                                         # "spring_meanTemp"))
# newpcdat <- cbind(pcdat, winter[, c("SiteID", "Year"), with=FALSE])
# newpcdat$Year <- as.character(newpcdat$Year)
# newpcdat <- newpcdat %>% 
#   mutate(zspring=scale_this(spring_meanTemp),
#          zwinter=scale_this(winter_meanTemp),
#          zfall=scale_this(fall_meanTemp),
#          zsummer=scale_this(prevsum_meanTemp))
winterdat$Year <- as.factor(as.character(winterdat$yearshift))
newdat2 <- left_join(newdat, winterdat, by = c("SiteID", "Year"))




#output pca results for table
pcout <- data.table(pc$rotation)
cols <- names(pcout)
pcout[,(cols) := round(.SD,2), .SDcols=cols]
pcout2 <- data.table(summary(pc)$importance)
pcout2[,(cols) := round(.SD,2), .SDcols=cols]

pcout <- rbind(pcout, pcout2)
pcout$Predictor <- c("Spring (t-1)",
                     "Summer (t-1)",
                     "Fall (t-1)",
                     "Winter",
                     "Spring",
                     "Standard deviation",
                     "Proportion of variance",
                     "Cumulative Proportion")
s <- setcolorder(pcout, names(pcout)[c(6, 1:5)])
write.csv(s, file = "pcatable.csv")
# 
# #pca shows promise for reducing variables to
# #spatial, temporal, and weird seasonal stuff
# pcdat <- newdat %>% 
#   dplyr::select(lat, meanGDD, winter_meanTemp,
#          fall_meanTemp, prevsum_meanTemp, prevspr_meanTemp) %>% 
#   distinct()
# pc <- prcomp(pcdat, center = TRUE, scale. = TRUE)
# pcdat <- cbind(pcdat, pc$x)
# newdat <- merge(newdat, pcdat, 
#                 by = c("lat", "meanGDD", "winter_meanTemp",
#                        "fall_meanTemp", "prevsum_meanTemp",
#                        'prevspr_meanTemp'))
newdat <- newdat2 %>% 
  group_by(species) %>%
  mutate(zDD = scale_this(log.lag.brood1),
         zlastbrood = scale_this(specscaleProp),
         DD2 = log.lag.brood1^2,
         zlat=scale_this(lat))
# test <- newdat %>% ungroup() %>% 
#   select(lat, zlastbrood, zsummer,zfall,zwinter,zspring)
#   t <- prcomp(test)
# covpca <- data.frame(t$x)
#   names(covpca) <- paste(names(covpca), "cov", sep="_")
# newdat1 <- cbind(as.data.frame(newdat),covpca)  
# testpc <- prcomp(newdat[,c("zlastbrood", "winter.onset", 
#                          "winter.end", "winter.length", 
#                          "winter.min", "winter.cummin", 
#                          "winter.snowtot", "winter.prectot")], 
#                center = TRUE, scale. = TRUE)
# test <- data.frame(testpc$x[,c(1:5)])
# names(test) <- c("npc1", "npc2", "npc3", "npc4", "npc5")
# newdat <- cbind(as.data.frame(newdat), test)

listcoef <- list()
listmod <- list()
listcoef2 <- list()
listmod2 <- list()
listplot <- list()
listimp <- list()
for (i in 1:length(modspec)){
  # SpeciesMod2 <- function(i){
    
    # scale_this <- function(x){
    #   (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    # }
    
  popdat <- newdat %>%
    filter(species == modspec[i]) # %>% 
    # mutate(zDD = scale_this(log.lag.brood1),
    #        zlastbrood = scale_this(specscaleProp),
    #        zlat = scale_this(lat),
    #        zfall = scale_this(fall_meanTemp),
    #        zwinter = scale_this(winter_meanTemp),
    #        zprevsum = scale_this(prevsum_meanTemp))
           # zgddleft = scale_this(gddleft),
           # zgddmismatch = scale_this(gddmismatch))
  # 
  # fit2 <- lm(lambda ~
  #              # log.lag.brood1 +
  #              sitescaleDD1 +
  #              PC1 +
  #                 PC2+
  #                 PC3+
  #                 PC4+
  #                 PC5
  #            ,
  #            data = popdat)

  fit2 <- lm(lambda ~
               # log.lag.brood1 +
               sitescaleDD1 +
               zlastbrood*
               (PC1 +
                  PC2+
                  PC3+
                  PC4+
                  PC5)
             ,
             data = popdat)
# 
#   fit2 <- lm(lambda ~
#               log.lag.brood1 +
#                # sitescaleDD1 +
#               zlastbrood *
#               (PC1 +
#                  PC2+
#                  PC3+
#                  PC4)
#             ,
#             data = popdat)
  # 
  # fit2 <- lmer(lambda ~ 
  #               sitescaleDD +
  #               specscaleProp + 
  #               I(specscaleProp^2) + 
  #               scalegddleft +
  #               I(scalegddleft^2) +
  #               (1|region9) +
  #               (1 + specscaleProp + 
  #                  I(specscaleProp^2) + 
  #                  scalegddmismatch +
  #                  I(scalegddmismatch^2)|Year),
  #             data = popdat)
  # 
  # outcoef <- tidy(fit)
  # outcoef$species <- modspec[i]
  # outmod <- glance(fit)
  # outmod$species <- modspec[i]
  outcoef2 <- tidy(fit2)
  outcoef2$species <- modspec[i]
  outmod2 <- glance(fit2)
  outmod2$species <- modspec[i]

  # 
  # listcoef[[i]] <- outcoef
  # listmod[[i]] <- outmod
  listcoef2[[i]] <- outcoef2
  listmod2[[i]] <- outmod2
  listplot[[i]] <- sjp.int(fit2,
                           show.ci = TRUE,
                           type = "eff",
                           title = modspec[i],
                           mdrt.values = "meansd",
                           geom.colors = c("blue", "dark green", "red"),
                           prnt.plot = FALSE)
  # listimp <- calc.relimp(fit2)
  
  # listimp[[i]] <- calc.relimp(fit2)
  # listimp[[i]] <- car::vif(fit2)
  # return(listimp)
  }
  
  
  library(parallel)
  # multicore
  system.time({
    cl <- makeCluster(8)
    clusterEvalQ(cl, {
      library(dplyr)
      library(relaimpo)
    })
    clusterExport(cl=cl,
                  varlist=c("newdat",
                            "modspec"))
    
    test <- parLapply(cl, 1:16, SpeciesMod2)
    stopCluster(cl)
  })

speccoef <- rbindlist(listcoef2)

cols <- names(speccoef)[2:5]
speccoef[,(cols) := round(.SD,3), .SDcols=cols]
speccoef %>% data.frame()
# specmod <- rbindlist(listmod)
specmod2 <- rbindlist(listmod2)
s <- speccoef %>% data.frame()
capture.output(s, file = "M2coef.final.txt")
s <- specmod2
capture.output(s, file = "M2glance.final.txt")
saveRDS(test, "M2relaimpoPCA.rds")
listimp <- readRDS("M2relaimpoPCA.rds")
# list of variable importance via relaimpo pkg
outdf <- data.frame()
for (i in 1:length(listimp)){
  species <- modspec[i]
  var.y <- listimp[[i]]@var.y
  r2 <- listimp[[i]]@R2
  lmg <- listimp[[i]]@lmg
  varnames <- names(lmg)
  avgcoef <- listimp[[i]]@ave.coeffs[, 10]
  outdf <- rbind(outdf, data.frame(species, var.y, r2, lmg, avgcoef, varnames))
}
row.names(outdf) <- NULL
s <- outdf  %>%  
  group_by(varnames) %>% 
  summarise(meanimp = mean(lmg / r2),
            sdimp = sd(lmg/r2),
            meanavgcoef = mean(avgcoef),
            sdavgcoef = sd(avgcoef)) %>% 
  arrange(-meanimp)
capture.output(s, file = "M2varimp.txt")

library(gridExtra)
#interaction plots from lm
#3 pdfs: 1 for each variable combo
int.grob <- list()
for (i in 1:16){
  int.grob[[i]] <- listplot[[i]]$plot.list[[1]]
}
n <- length(int.grob)
nCol <- floor(sqrt(n))
pdf("M2interaction1pca.pdf", height=15, width=20)
do.call("grid.arrange", c(int.grob, ncol=nCol))
dev.off()


#some interactions significant when all species together
#lat/fall temp interact with extra brood


# defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e8)

# trouble getting convergence with varying slopes by species

# species lastbrood with regional centering
specdat <- newdat %>% 
  group_by(species, region9) %>% 
  mutate(zbroodreg = scale_this(specscaleProp),
         varbrood = specscaleProp - mean(specscaleProp),
         meanbrood = mean(specscaleProp),
         varbrood2 = zlastbrood - mean(zlastbrood),
         meanbrood2 = mean(zlastbrood))
         # zPC1reg = scale_this(PC1),
         # zPC2reg = scale_this(PC2),
         # zPC3reg = scale_this(PC3),
         # zPC4reg = scale_this(PC4),
         # PC1var = PC1 - mean(PC1),
         # PC1mean = mean(PC1),
         # PC2var = PC2 - mean(PC2),
         # PC2mean = mean(PC2),
         # PC3var = PC3 - mean(PC3),
         # PC3mean = mean(PC3),
         # PC4var = PC4 - mean(PC4),
         # PC4mean = mean(PC4)) 
# %>% 
#   ungroup() %>% 
#   group_by(species, Year) %>% 
#   mutate(zbroodyr = scale_this(specscaleProp),
#          zPC1yr = scale_this(PC1),
#          zPC2yr = scale_this(PC2),
#          zPC3yr = scale_this(PC3),
#          zPC4yr = scale_this(PC4))

bsdat <- specdat %>% 
  filter(species == "Spicebush Swallowtail") #%>%
  # mutate(PC1 = scale_this(PC1),
         # PC2 = scale_this(PC2),
         # PC3 = scale_this(PC3),
         # PC4 = scale_this(PC4))

t <- ggplot(bsdat, aes(x=PC1var, y = varbrood, fill = region9, color = region9)) +
  geom_point()
t

#problem with varying slopes by region, 
#patterns don't fit expectations, no clear gradient or even
#spatial correlation between regional effects

# fit0 <- lmer(lambda ~
#                sitescaleDD1 +
#                +specscaleProp
#              *(PC1
#                +PC2
#                +PC3
#                +PC4)
#              +(0
#                +specscaleProp
#                +PC1
#                +PC2
#                +PC3
#                +PC4
#                +specscaleProp:PC1
#                +specscaleProp:PC2
#                +specscaleProp:PC3
#                +specscaleProp:PC4
#                |region9)
#           ,
#              data = bsdat, REML = TRUE,
#              control = lmerControl(optCtrl = list(maxfun = 1e8)))

fit1 <- lmer(lambda ~
               sitescaleDD1 +
               +zlastbrood
             *(PC1
               +PC2
               +PC3
               +PC4
               )
             +(0
               +zlastbrood
               # *PC1
               |region9)
             ,
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))




fit2 <- lmer(lambda ~
                sitescaleDD1
              +varbrood2*meanbrood2
              +site.var.pc1*site.mean.pc1
              +site.var.pc2*site.mean.pc2
              +site.var.pc3*site.mean.pc3
              +site.var.pc4*site.mean.pc4
              +varbrood2:site.var.pc1
             +varbrood2:site.var.pc2
             +varbrood2:site.var.pc3
             +varbrood2:site.var.pc4
             
              +
                (0
                 +varbrood2
                 # +site.var.pc1
                 # +varbrood2:site.var.pc1
                 |region9)
       , 
              data = bsdat, REML = TRUE,
              control = lmerControl(optCtrl = list(maxfun = 1e8)))





fit3 <- lmer(lambda ~
               sitescaleDD1
               +zlastbrood
             +site.var.pc1*site.mean.pc1
             +site.var.pc2*site.mean.pc2
             +site.var.pc3*site.mean.pc3
             +site.var.pc4*site.mean.pc4
             +zlastbrood:site.var.pc1
             +zlastbrood:site.var.pc2
             +zlastbrood:site.var.pc3
             +zlastbrood:site.var.pc4
              +
               (0 
                +zlastbrood
                |region9), 
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))

# bsdat$pred <- predict(fit2)
# ggplot(data = bsdat, aes(x = site.mean.pc1, y = site.var.pc1,
#                          fill = pred, color = pred)) + stat_density()

AIC(fit1, fit2, fit3)
rsquared.glmm(fit1)
rsquared.glmm(fit2)
rsquared.glmm(fit3)


fit3 <- lmer(lambda ~
               sitescaleDD1
             +varbrood*
               (meanbrood
                +PC1
                +PC2
                +PC3
                +PC4
               )             
             # + (1|Year)
              + (1 
                 +zlastbrood
                |region9), 
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))


fit6 <- lmer(lambda ~
               sitescaleDD1
             +varbrood*(meanbrood
             +site.var.pc1*site.mean.pc1
             +site.var.pc2*site.mean.pc2
             +site.var.pc3*site.mean.pc3
             +site.var.pc4*site.mean.pc4)
             +
               (1 
                |region9), 
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))


mod0 <- lmer(lambda ~ 
               sitescaleDD1 
             + (1|Year)     ,
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))
# control=lmerControl(optimizer="nloptwrap2"))
mod1 <- lmer(lambda ~ 
                       sitescaleDD1 +
             zlastbrood *
               (PC1 +
                  PC2 +
                  PC3 +
                  PC4
               ) 
                     + (1|Year)     ,
                     data = bsdat, REML = TRUE,
                     control = lmerControl(optCtrl = list(maxfun = 1e8)))

mod2 <- lmer(lambda ~ 
               sitescaleDD1 +
               zlastbrood *
               (PC1 +
                  PC2 +
                  PC3 +
                  PC4
               ) 
             + (1|region9)     ,
             data = bsdat, REML = TRUE,
             control = lmerControl(optCtrl = list(maxfun = 1e8)))


bsmod2a <- lmer(lambda ~ 
       sitescaleDD1 +
       zlastbrood *
       (zPC1reg +
          zPC2reg +
          zPC3reg +
          zPC4reg
       ) 
       + (1|Year)     ,
     data = bsdat, REML = TRUE,
     control = lmerControl(optCtrl = list(maxfun = 1e8)))
     # control=lmerControl(optimizer="nloptwrap2"))

bsmod2 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 zlastbrood *
                 (PC1 +
                    PC2 +
                    PC3 +
                    PC4
                 ) 
                  + (0 
                     + zlastbrood 
                     + PC1
                     +PC2
                     +PC3
                     +PC4
                     +zlastbrood:PC1
                     +zlastbrood:PC2
                     +zlastbrood:PC3
                     +zlastbrood:PC4
                     |species)
               + (1|Year)     ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))

bsmod2 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood2*
                 (PC1var +
                    PC2var +
                    PC3var +
                    PC1mean +
                    PC2mean +
                    PC3mean)+
                 meanbrood2 +
               + (1|Year),
               data = bsdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))

bsmod2 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 # varbrood+
                 # PC1var +
                    # PC2var +
                    PC3var+
                 # meanbrood +
                 # PC1mean+
                    # PC2mean +
                    PC3mean+
                 + (1|Year),
               data = bsdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))

# mod <- gam(lambda~s(zlastbrood) +
#              s(PC1) + ti(zlastbrood, PC1), data = bsdat)

#trying to do each PC separately, to include lots of interactions
#trouble fitting some, rescale suggested, 3-way interaction tough to grok

fitpc1 <- lmer(lambda ~ 
                   sitescaleDD1 +
                   varbrood*
                   meanbrood
                 +PC1var*PC1mean
                +varbrood:PC1var
               +varbrood:PC1mean
               +varbrood:PC1var:PC1mean
                 +(1
                   +sitescaleDD1
                   +varbrood
                   +meanbrood
                   +varbrood:meanbrood
                   +PC1var
                   +PC1mean
                   +PC1var:PC1mean
                   +varbrood:PC1var
                   +varbrood:PC1mean
                   +varbrood:PC1var:PC1mean
                   
                   |species)
                 +
                   (1|Year)
                 ,
                 data = specdat, REML = TRUE,
                 # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                 control=lmerControl(optimizer="nloptwrap2"))

fitPC2 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood*
                 meanbrood
               +PC2var*PC2mean
               +varbrood:PC2var
               +varbrood:PC2mean
               +varbrood:PC2var:PC2mean
               
               +(1
                 +sitescaleDD1
                 +varbrood
                 +meanbrood
                 +varbrood:meanbrood
                 +PC2var
                 +PC2mean
                 +PC2var:PC2mean
                 +varbrood:PC2var
                 +varbrood:PC2mean
                 +varbrood:PC2var:PC2mean
                 
                 |species)
               +
                 (1|Year)
               ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))


fitPC3 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood*
                 meanbrood
               +PC3var*PC3mean
               +varbrood:PC3var
               +varbrood:PC3mean
               +varbrood:PC3var:PC3mean
               
               +(1
                 +sitescaleDD1
                 +varbrood
                 +meanbrood
                 +varbrood:meanbrood
                 +PC3var
                 +PC3mean
                 +PC3var:PC3mean
                 +varbrood:PC3var
                 +varbrood:PC3mean
                 +varbrood:PC3var:PC3mean
                 
                 |species)
               +
                 (1|Year)
               ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))


fitPC4 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood*
                 meanbrood
               +PC4var*PC4mean
               +varbrood:PC4var
               +varbrood:PC4mean
               +varbrood:PC4var:PC4mean
               
               +(1
                 +sitescaleDD1
                 +varbrood
                 +meanbrood
                 +varbrood:meanbrood
                 +PC4var
                 +PC4mean
                 +PC4var:PC4mean
                 +varbrood:PC4var
                 +varbrood:PC4mean
                 +varbrood:PC4var:PC4mean
                 
                 |species)
               +
                 (1|Year)
               ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))

all0 <- lmer(lambda ~
               sitescaleDD1
             +varbrood*meanbrood
             +site.var.pc1*site.mean.pc1
             +site.var.pc2*site.mean.pc2
             +site.var.pc3*site.mean.pc3
             +site.var.pc4*site.mean.pc4
             +
               (1 
                +sitescaleDD1
                +varbrood*meanbrood
                # +site.var.pc1*site.mean.pc1
                # +site.var.pc2*site.mean.pc2
                # +site.var.pc3*site.mean.pc3
                # +site.var.pc4*site.mean.pc4
                |species)
             +(1|SiteYear),
              
             data = specdat, REML = TRUE,
             control=lmerControl(optimizer="nloptwrap2"))


all1d <- lmer(lambda ~
       sitescaleDD1
     +varbrood2*meanbrood2
     +site.var.pc1*site.mean.pc1
     +varbrood2:site.var.pc1
     +site.var.pc2*site.mean.pc2
     +varbrood2:site.var.pc2
     +(0
       +sitescaleDD1
       +varbrood2
       # +meanbrood2
       +varbrood2:meanbrood2
       +site.var.pc1
       +site.var.pc2
       # +site.mean.pc1
       # +site.var.pc1:site.mean.pc1
       +varbrood2:site.var.pc1 
       +varbrood2:site.var.pc2
       |species)
     +(1|SiteYear)
     ,
     data = specdat, REML = TRUE,
     control=lmerControl(optimizer="nloptwrap2"))

all3ae <- lmer(lambda ~
               sitescaleDD1
             +(  site.var.pc1
                 +site.var.pc2
                 +site.var.pc3
                 +site.var.pc4
               )
             *varbrood2
             +
               (0
                +sitescaleDD1
                +site.var.pc1
                +site.var.pc2
                +site.var.pc3
                +site.var.pc4
                +varbrood2
                +varbrood2:site.var.pc1
                +varbrood2:site.var.pc2
                +varbrood2:site.var.pc3
                +varbrood2:site.var.pc4
                |species) 
             +(1|SiteYear)
             ,
             data = specdat, REML = TRUE,
             control=lmerControl(optimizer="nloptwrap2"))

AIC(all3ac, all3ad, all3ae)

all3b <- lmer(lambda ~
               sitescaleDD1
             +(site.var.pc1+site.mean.pc1
               +site.var.pc2+site.mean.pc2
               +site.var.pc3+site.mean.pc3
               +site.var.pc4+site.mean.pc4)
             *zlastbrood
             +
               (1
                # +sitescaleDD1
                # +site.var.pc1+site.mean.pc1
                # +site.var.pc2+site.mean.pc2
                # +site.var.pc3+site.mean.pc3
                # +site.var.pc4+site.mean.pc4
                # +zlastbrood
                # +zlastbrood:PC1
                # +zlastbrood:PC2
                # +zlastbrood:PC3
                # +zlastbrood:PC4
                |species) 
             +(1|SiteYear),
             data = specdat, REML = TRUE,
             control=lmerControl(optimizer="nloptwrap2"))




fitalla4 <- lmer(lambda ~ 
                   sitescaleDD1 +
                   varbrood+
                     meanbrood
                   +PC1var+PC1mean
                 +PC2var+PC2mean
                 +PC3var+PC3mean
                 +PC4var+PC4mean
                 
                    +
                   (0
                    +varbrood
                    +PC1var
                    +PC2var
                    +PC3var
                    +PC4var
                    |region9) 
                 +(1
                   +sitescaleDD1
                   +varbrood
                   +PC1var
                   +PC2var
                   +PC3var
                   +PC4var
                   +meanbrood
                     +PC1mean
                     +PC2mean
                   +PC3mean
                   +PC4mean|species)
                 +
                 (1|Year)
                 ,
                 data = specdat, REML = TRUE,
                 # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                 control=lmerControl(optimizer="nloptwrap2"))

fitall2 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood2 +
                  meanbrood2
                   +PC1var*PC1mean
                   +PC2var*PC2mean
                   +PC3var*PC3mean
                   +PC4var*PC4mean
               +
                 (1
                  +sitescaleDD1
                  +varbrood2
                  # +meanbrood
                  +PC1var
                  +PC2var
                  +PC3var
                  +PC4var
                  # +varbrood:PC1var
                  # +varbrood:PC2var
                  # +varbrood:PC3var
                  |species) 
               +
                 (1|Year)
               ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))


fitall3 <- lmer(lambda ~ 
                 sitescaleDD1 +
                 varbrood2*
                   meanbrood2*(
                   +PC1
                   # +PC2
                   # +PC3
                   # +PC4
                   )
               +
                 (0
                  +sitescaleDD1
                  +varbrood2
                  +meanbrood2
                  +varbrood2:meanbrood2
                  +PC1
                  # +PC2
                  # +PC3
                  # +PC4
                  + varbrood2:PC1
                  # + varbrood2:PC2
                  # + varbrood2:PC3
                  # + varbrood2:PC4
                  + meanbrood2:PC1
                  + varbrood2:meanbrood2:PC1
                  |species) 
               +
                 (1|SiteYear)
               ,
               data = specdat, REML = TRUE,
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))
               control=lmerControl(optimizer="nloptwrap2"))


fitall1a <- lmer(lambda ~ 
                  sitescaleDD1 +
                  zbroodreg *
                  (PC1var +
                     PC2var+
                     PC3var
                  ) +
                  (1
                   +sitescaleDD1
                   +zbroodreg
                   +PC1var
                   +PC2var
                   +PC3var
                   +zbroodreg:PC1var
                   +zbroodreg:PC2var
                   +zbroodreg:PC3var
                   |species) + 
                  (1|Year)
                ,
                data = specdat, REML = TRUE,
                # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                control=lmerControl(optimizer="nloptwrap2"))



fitb <- lmer(lambda ~
                   sitescaleDD1 +
                   varbrood2 *
                   (PC1 +
                      PC2+
                      PC3 +
                      PC4
                   ) +
                  meanbrood2+
                   (1
                    +sitescaleDD1
                    +varbrood2
                    +meanbrood2
                    +PC1
                    +PC2
                    +PC3
                    +PC4
                    +varbrood2:PC1
                    +varbrood2:PC2
                    +varbrood2:PC3
                    +varbrood2:PC4
                    |species)
                 ,
                 data = specdat, REML = TRUE,
                 # control = lmerControl(optCtrl = list(maxfun = 1e8)))
control=lmerControl(optimizer="nloptwrap2"))

AIC(fita, fitb)

specdat$SiteYrID <- paste(specdat$SiteID, specdat$Year, sep = "_")
newdat$SiteYrID <- paste(newdat$SiteID, newdat$Year, sep = "_")
newdat$RegYrID <- paste(newdat$region9, newdat$Year, sep = "_")

specdat <- newdat %>% 
  group_by(species, region9) %>% 
  mutate(
         npc1var = npc1 - mean(npc1),
         npc1mean = mean(npc1),
         npc2var = npc2 - mean(npc2),
         npc2mean = mean(npc2),
         npc3var = npc3 - mean(npc3),
         npc3mean = mean(npc3),
         npc4var = npc4 - mean(npc4),
         npc4mean = mean(npc4),
         npc5var = npc5 - mean(npc5),
         npc5mean = mean(npc5) )



newfit1 <- lmer(lambda ~
                  sitescaleDD1
                +(PC1
                +PC2
                +PC3
                +PC4)
                *varbrood2
                +
                  (1
                   # +sitescaleDD1
                   +PC1
                   +PC2
                   +PC3
                   +PC4
                   +varbrood2
                   # +zlastbrood:PC1
                   # +zlastbrood:PC2
                   # +zlastbrood:PC3
                   # +zlastbrood:PC4
                   |species/region9),
                data = specdat, REML = TRUE,
                control=lmerControl(optimizer="nloptwrap2"))

newfit1 <- lmer(lambda ~
                 sitescaleDD1
               +npc1
               +npc2
               +npc3
               +npc4
               +npc5
               +
                 (1
                  # +sitescaleDD1
                  # +npc1
                  # +npc2
                  # +npc3
                  # +npc4
                  # +npc5
                  |species)
               , 
               data = newdat, REML = TRUE,
               control=lmerControl(optimizer="nloptwrap2"))
                # control=lmerControl(optimizer="nloptwrap2",
                #     optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
               # control = lmerControl(optCtrl = list(maxfun = 1e8)))

bsdat <- specdat %>% 
  filter(species == "Spicebush Swallowtail")

newfit5 <- lmer(lambda ~
                  sitescaleDD1
                +npc1var
                +npc2var
                +npc3var
                +npc4var
                +npc5var
                +npc1mean
                +npc2mean
                +npc3mean
                +npc4mean
                +npc5mean
                +npc1var:npc1mean
                +npc2var:npc2mean
                +npc3var:npc3mean
                +npc4var:npc4mean
                +npc5var:npc5mean
                +
                  (1
                   # +npc1var
                   # +npc2var
                   # +npc3var
                   # +npc4var
                   # +npc5var
                   |species) + (1|Year)
                , 
                data = specdat, REML = TRUE,
                control=lmerControl(optimizer="nloptwrap2"))
# control=lmerControl(optimizer="nloptwrap2",
#     optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
# control = lmerControl(optCtrl = list(maxfun = 1e8)))


newfit5b <- lmer(lambda ~
                  sitescaleDD1
                +npc1var
                  +npc2var
                  +npc3var
                  +npc4var
                  +npc5var
                +npc1mean
                +npc2mean
                +npc3mean
                +npc4mean
                +npc5mean
                +
                  (1
                   +sitescaleDD1
                   +npc1var
                   +npc2var
                   +npc3var
                   +npc4var
                   +npc5var
                   +npc1mean
                   +npc2mean
                   +npc3mean
                   +npc4mean
                   +npc5mean
                   |species)
                
                +(1
                  |SiteYrID)
                , 
                data = specdat, REML = TRUE,
                control=lmerControl(optimizer="nloptwrap2"))
# control=lmerControl(optimizer="nloptwrap2",
#     optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
# control = lmerControl(optCtrl = list(maxfun = 1e8)))


specdat$species <- as.factor(specdat$species)
m <- gam(lambda ~
           s(sitescaleDD1)
         +s(varbrood2)
         +s(meanbrood2)
         +s(PC1var, by = species)
         +s(PC2var)
         +s(PC3var)
         +s(PC4var), data = specdat)



fitb <- lmer(lambda ~
                   sitescaleDD1 +
                   varbrood2 *
                   (PC1 +
                      PC2+
                      PC3 +
                      PC4
                   ) +
                   meanbrood2+
                   (1
                    # +sitescaleDD1
                    +varbrood2
                    +meanbrood2
                    # +PC1
                    # +PC2
                    # +PC3
                    # +PC4
                    # +varbrood2:PC1
                    # +varbrood2:PC2
                    # +varbrood2:PC3
                    # +varbrood2:PC4
                    |species)
             +(1|SiteYrID)
                 ,
                 data = specdat, REML = TRUE,
                 # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                 control=lmerControl(optimizer="nloptwrap2"))

AIC(fita, fitb)

# saveRDS(g0.bobyqa2, "M2lmer.fit.pca.final.rds")
saveRDS(fitall2a, "M2lmer.final.rds")

fitall2a <- readRDS("M2lmer.final.rds")
s <- summary(fitall2a)
capture.output(s, file = "M2lmer.final.txt")

sjt.lmer(fitall2a, p.kr = FALSE, file = "M2table.html",
         show.ci = FALSE, show.se = TRUE, digits.est = 3,
         digits.se = 3)

sjt.lmer(fitall2a, 
  pred.labels = c("Previous year density", "Size of extra brood", 
                  "PC1 (Warmer in all seasons)",
                  "PC2 (Warm winter)", "PC3 (Cool fall, warm spring)", 
                  "PC4 (Cool prev spring, warm summer)", "PC5 (Cool spring/fall, warm prev spring)",
                  "Extra brood x PC1", "Extra brood x PC2", "Extra brood x PC3",
                  "Extra brood x PC4", "Extra brood x PC5"),
  depvar.labels = "Response: Annual growth rate",
  p.kr = FALSE,
  show.ci = FALSE, show.se = TRUE,
  show.icc = FALSE, show.re.var = FALSE,
  string.se = "std. error", cell.spacing = .05,
  digits.est = 3, digits.se = 3, file = "m2fit.html")
  
null <- gls(lambda ~
                   sitescaleDD1  +
                   zlastbrood *
                   (PC1 +
                      PC2+
                      PC3+
                      PC4 +
                      PC5),
            data = newdat,
            method = "REML")

#try removing pearley eye
newdat <- newdat %>% 
  filter(species != "Northern Pearly-Eye")

  fitall2a <- lmer(lambda ~ 
                                sitescaleDD1 +
                                zlastbrood *
                                (PC1 +
                                   PC2+
                                   PC3+
                                   PC4 + 
                                   PC5) +
                                (1
                                 +sitescaleDD1
                                 +zlastbrood
                                 +PC1
                                 +PC2
                                 +PC3
                                 +PC4
                                 +PC5
                                 +zlastbrood:PC1
                                 +zlastbrood:PC2
                                 +zlastbrood:PC3
                                 +zlastbrood:PC4
                                 +zlastbrood:PC5
                                 |species)
                              ,
                              data = newdat, REML = TRUE,
                              # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                   control=lmerControl(optimizer="nloptwrap2"))
  
  AIC(fitall2a, fitall2b)
  # AIC(null, fitall2a)
  
fitall2b <- lmer(lambda ~ 
                   sitescaleDD1 +
                  zlastbrood *
                  (PC1 +
                     PC2+
                     PC3+
                     PC4+
                     PC5
                     ) +
                  (1
                   +sitescaleDD1
                  +zlastbrood
                   +PC1
                   +PC2
                   +PC3
                   +PC4
                   +PC5
                   +zlastbrood:PC1
                   +zlastbrood:PC2
                   +zlastbrood:PC3
                   +zlastbrood:PC4
                   +zlastbrood:PC5
                   |species)
                ,
                data = newdat, REML = TRUE,
                # control = lmerControl(optCtrl = list(maxfun = 1e8)))
                control=lmerControl(optimizer="nloptwrap2"))

AIC(fitall2a, fitall2b)



#convergence issues

library(nloptr)
## from https://github.com/lme4/lme4/issues/98:
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e7)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
  for (n in names(defaultControl)) 
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}
g0.bobyqa2 <- update(fitall2a,control=lmerControl(optimizer="nloptwrap2"))
g0.NM2 <- update(fitall2,control=lmerControl(optimizer=nloptwrap2,
                                               optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))

# plot part 2 lmer results: species variation
preddat <- newdat %>% 
  group_by(species) %>% 
  summarise_at(vars(PC1,PC2,PC3,PC4,PC5),
               funs(mean, sd))


preddata <- expand.grid(species = unique(newdat$species),
                      sitescaleDD1 = 0,
                      zlastbrood = c(-1, 1),
                      PC1 = 0,
                      PC2 = 0,
                      PC3 = 0,
                      PC4 = 0,
                      PC5 = seq(-1.15,2.25,0.05))
preddata <- merge(preddat,preddata)
# only use PC data for each species within 1 SD of its mean
preddata <- preddata %>%
  rowwise() %>% 
  mutate(PCmin = PC5_mean - PC5_sd,
         PCmax = PC5_mean + PC5_sd) %>% 
  filter(PC5 >= PCmin && PC5 <= PCmax)
preddata$lambda <- 0


fm1 <- fitall2a
mm <- model.matrix(terms(fm1), preddata)
preddata$lambda <- predict(fm1, preddata, re.form=NA)
preddata$lambda2 <- predict(fm1, preddata, re.form=NULL)
preddata$zlastbrood <- as.factor(as.character(preddata$zlastbrood))
levels(preddata$zlastbrood) <- c("Smaller extra brood", "Larger extra brood")
alldat <- distinct(preddata[, c("zlastbrood", "PC5", "lambda")])

specdat <- preddata %>% 
  group_by(species, zlastbrood) %>% 
  summarise(PC5 = max(PCmax))

library(ggrepel)
plt <- ggplot(data=preddata, aes(x=PC5, y=lambda2, 
                                 group = interaction(species, zlastbrood))) +
  geom_line(alpha=.3, size = .8, color = "black") + 
  facet_grid(~zlastbrood) +
  # geom_text_repel(
  #   data = specdat,
  #   aes(label = species),
  #   size = 6) +
  theme(legend.position = "none") +
  labs(x = "PC5: Warm prev. spring, cool winter/spring",
       y = "Annual growth rate") +
  geom_line(data=alldat, aes(x=PC5, y = lambda, group = zlastbrood), color = "black", size = 1.5) +
  theme_bw()

pdf(paste("M2plotPC5.pdf", sep = ""), width = 4, height = 2)
print(plt)
dev.off()



# predict lambda under different scenarios
library(merTools)
mod <- fitall2a
plotREsim(REsim(mod))

exam <- draw(mod, type ="average", 
             varList=list("species"="Black Swallowtail"))
exam <- wiggle(exam, var = "zlastbrood", 
                   values = c("0", "1"))

exam$yhat <- predict(mod, newdata = exam)

ggplot(example3, aes(x = service, y = yhat)) + geom_line(aes(group = 1)) + 
  theme_bw() + ylim(c(1, 5)) +
  geom_hline(yintercept = mean(InstEval$y), linetype = 2) + 
  geom_hline(yintercept = mean(InstEval$y) + sd(InstEval$y), linetype = 3) + 
  geom_hline(yintercept = mean(InstEval$y) - sd(InstEval$y), linetype = 3)



coefpca <- coef(fitall2a)$species
coefpca$species <- rownames(coefpca)
df <- coefpca[,c(3:13)]
m2pc <- prcomp(df, center = TRUE, scale. = TRUE)
m2pcdf <- as.data.frame(m2pc$x)
names(m2pcdf) <- paste(names(m2pcdf), "coef", sep = "_")
coefpca <- cbind(coefpca, m2pcdf)

a <- ggplot(data=coefpca, aes(x=PC1_coef,y=PC2_coef,label=species))+
  geom_text()
a

newdat$species <- as.factor(newdat$species)
# not much added with gam
gmod <- gam(lambda ~ s(sitescaleDD1) + 
               s(zlastbrood) + 
               s(scalefall) +
               s(scalewinter) +
              s(Year, bs = "re")+
              s(species,bs="re"),
             data = newdat)


specdat <- newdat %>% 
  filter(species=="Spicebush Swallowtail") %>% 
  group_by(region9) %>% 
  mutate(vlastbrood = lag.lastprop - mean(lag.lastprop),
         mlastbrood = mean(lag.lastprop),
         vfall = zfall - mean(zfall),
         vwinter = zwinter - mean(zwinter),
         mfall = mean(zfall),
         mwinter = mean(zwinter))
mod2 <- lmer(lambda ~ sitescaleDD1 +
             vlastbrood *
             (vfall+
                vwinter) +
               mlastbrood +
               mfall +
               mwinter +
             (1 +
                vfall+
                   vwinter
                |Year),
           data = specdat, REML = TRUE,
           # control = lmerControl(optCtrl = list(maxfun = 1e8)))
           control=lmerControl(optimizer="nloptwrap2"))
             

#these give neat visualizations of nonlinear effects of
#latitude, temperature, ordinal date
#possible good supplement to linear models.

library(randomForest)
library(forestFloor)

rflist <- list()
for (i in 1:length(modspec)){
  spec <- modspec[i]
  print(spec)
  rfdat <- newdat[which(newdat$species == spec), ]
  rfdat <- rfdat[, names(newdat)[c(26:28, 30, 39, 40)]]
  rfdat <- rfdat[complete.cases(rfdat), ]
  
  y <- rfdat$log.brood1
  x <- rfdat[, -3]
  
  mod <- randomForest(x = x, y = y, 
                      importance = TRUE, 
                      ntree = 1001, 
                      mtry = 2,
                      keep.inbag = TRUE,
                      keep.forest = TRUE)
  rflist[[i]] <- mod
  rflist[[i]]$species <- spec
}

#compute forestFloor object, often only 5-10% time of growing forest
ff = forestFloor(
  rf.fit = mod,       # mandatory
  X = x,              # mandatory
  calc_np = FALSE,    # TRUE or FALSE both works, makes no difference
  binary_reg = FALSE  # takes no effect here when rfo$type="regression"
)

Col = fcol(ff,cols=2, outlier.lim = 2)

#plot partial functions of most important variables first
plot(ff,                       # forestFloor object
     # plot_seq = ,           # optional sequence of features to plot
     orderByImportance=TRUE,
     plot_GOF = TRUE,
     col = Col
)

library(edarf)




library(lme4)
mod3 <- lmer(lambda ~ 
              (sitescaleDD1 +
               sitescaleDDall +
                scalegddmismatch +
                scalegddleft +
                scalewinter +
                scalefall) *
               specscaleProp +
               I(specscaleProp^2) +
              (1 + specscaleProp + I(specscaleProp^2) |species) +
               (1|Year),
            data = newdat)

mod <- lmer(lambda ~ 
              
              sitescaleDDall +
              specscaleProp + 
                 scalegddleft +
                 scalewinter +
              (1|Year) + 
              (1 + sitescaleDDall +
                 specscaleProp + 
                    scalegddleft +
                    scalewinter|species),
            data = newdat)

mod <- lmer(lambda ~ 
              sitescaleDD +
              specscaleProp * 
              scalegddleft +
              (1|Year),
            data = popdat)

summary(mod)
sjp.int(fit, type = "eff")  
r.squared.merMod(mod3)

# 
# mod <- lmer(lambda ~ sitescaleDD + 
#               lastprop * 
#               lat + 
#               (1|Year) + (1|species),
#             data = newdat)



broodplot <- ggplot(newdat, aes(x = lon, y = lat, color = scaleProp)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Year)
broodplot


mod1 <- lmer(lambda ~ sitescaleDD +scalelat * ( 
  scaleprop + I(scaleprop^2)) +
    (1|Year),
  data = newdat)

mod2 <- lmer(lambda ~ sitescaleDD + 
               lastprop + I(lastprop^2) +
               (1|Year),
             data = newdat)

library(effects)
eff <- allEffects(fitall)
plot(eff)


library(sjPlot)
# sjp.setTheme(theme = theme_minimal())

sjt.lmer(mod1, file = "brood.html")


sjp.lmer(fitall,facet.grid = FALSE,
         sort.coef = "sort.all")
sjp.lmer(fitall, type = "resp")

sjp.lmer(fitall, type = "poly", poly.term = "scalewinter")
sjp.int(fitall, type = "cond")  


library(effects)
eff <- allEffects(fit1)
eff <- Effect("ztemp", fit1)
plot(eff)

# not much added with gam
gmod <- gamm(lambda ~ s(scaleDD) + 
               s(scaleProp) + 
               s(lat), random = list(Year=~1),
             data = newdat)


