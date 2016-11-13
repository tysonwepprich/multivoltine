rm(list=ls())
library(mclust)
library(plyr)
library(dplyr)
library(mgcv)
library(lubridate)
library(tidyr)
library(stringr)
library(ggplot2)
library(geosphere)
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

library(data.table)
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
gdd <- gdd[gdd$year >= 1995, ]
gdd <- as.data.frame(gdd)
gdd$date <- strptime(paste(gdd$year, gdd$yday, sep = " "), format = "%Y %j")
gdd$SiteDate <- as.Date(gdd$date)
gdd$date <- NULL

gdd <- gdd %>%
  mutate(SiteID = formatC(as.numeric(site), width = 3, format = "d", flag = "0"))
gdd <- merge(gdd, site_geo, by = "SiteID")
gdd$SiteYear <- paste(gdd$SiteID, gdd$year, sep = "_")
gdd_tomerge <- gdd %>% select(SiteYear, yday, cumdegday) %>%
  dplyr::rename(Ordinal = yday)
rm(gdd)

# i <- 67
fs <- list.files("gamGDD")
#only try a few species first for mixture model brood separation
mvspec <- c(10, 24, 28, 43, 30, 47, 49, 56, 59, 63, 67, 80, 86, 77, 83, 71)
mvbroodmin <- c(2, 1, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 3, 2, 2, 2)
mvbroodmax <- c(3, 2, 3, 4, 2, 3, 2, 2, 2, 3, 3, 3, 4, 3, 3, 3)

# worried about uncertainty and bad clustering when low populations at sites
counts <- list()
for(i in 1:16){
  fs <- list.files("gamGDD")
  gamlist <- readRDS(paste("gamGDD/", fs[mvspec[i]], sep = ""))
  # gammod <- gamlist$mod
  # pred <- gamlist$preds
  temp <- gamlist$counts
  rm(gamlist)
  temp$filenum <- i
  counts[[i]] <- temp
}

allcounts <- rbindlist(counts)
allcounts$species <- fs[mvspec[allcounts$filenum]]
allcounts$species <- gsub('gam', '', allcounts$species)
allcounts$species <- gsub('.rds', '', allcounts$species)
siteyrcounts <- allcounts %>% 
  group_by(species, SiteYear) %>% 
  summarise(tots = sum(Total))
regsites <- allcounts %>% 
  group_by(Reg9Year) %>% 
  summarise(numsites = length(unique(SiteID)))


# for (i in 1:length(fs)){
BroodMixMod <- function(i){
  fs <- list.files("gamGDD")
  gamlist <- readRDS(paste("gamGDD/", fs[mvspec[i]], sep = ""))
  gammod <- gamlist$mod
  pred <- gamlist$preds
  datGAM <- gamlist$counts
  rm(gamlist)
  
  newpred <- pred %>%
    data.frame() %>%
    dplyr::select(SiteID, SiteYear, Reg9Year) %>%
    distinct() 
  
  #GDD prediction by arbitrary binwidth
  # allregs <- expand.grid(unique(newpred$Reg9Year), 
  #                        seq(0, max(pred$cumdegday), 5))
  # names(allregs) <- c("Reg9Year", "cumdegday")
  # pred <- full_join(newpred, allregs, by = "Reg9Year", all.x = TRUE, all.y = TRUE)
  # pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response")) 
  #                                     # exclude = c("s(SiteID)"))
  # pred <- as.data.frame(pred)
  
  #GDD prediction, but by ordinal day (butterfly day index)
  allregs <- expand.grid(unique(newpred$Reg9Year),
                         seq(91, 305, 1))
  names(allregs) <- c("Reg9Year", "Ordinal")
  pred <- full_join(newpred, allregs, by = "Reg9Year", all.x = TRUE, all.y = TRUE)
  pred <- left_join(pred, gdd_tomerge, by = c("SiteYear", "Ordinal"))
  
  
  pred$GAM.pred.reg9yr <- as.vector(predict.gam(gammod, pred, type = "response"))
  pred$GAM.pred.reg9yr2 <- as.vector(predict.gam(gammod, pred, type = "response",
                                                exclude = c("s(SiteID)")  ))
  pred <- as.data.frame(pred)
  
  
  nsim <- 100
  # outlist <- list()
  outdf <- data.frame()
  for(ns in 1:nsim){
    for (j in 1:length(unique(pred$SiteYear))){
      
      index <- unique(pred$SiteYear)[j]
      clustdat <- pred %>%
        filter(SiteYear == index) %>%
        # select(Ordinal, GAM.pred, lat, SiteID) %>%
        dplyr::select(Ordinal, cumdegday, GAM.pred.reg9yr) %>%
        rowwise() %>%
        mutate(GAMpois = rpois(1, GAM.pred.reg9yr[1])) %>%
        ungroup()
      # 
      # dumbFUNC <- function(x) sum(rpois(100, x))
      # clustdat$GAMpois <- sapply(clustdat$GAM.pred.reg9yr, FUN = dumbFUNC)
      # 
      if(sum(clustdat$GAMpois) < 10){
        next
      }
      
      dd <- rep(clustdat$cumdegday, clustdat$GAMpois)
      # daycount <- rep(clustdat$Ordinal, clustdat$GAMpois)
      # ddcount <- rep(clustdat$cumdegday, clustdat$GAMpois)
      
      # dat <- data.frame(cbind(yd, lat, dd))
      # dat <- data.frame(cbind(lat, dd))
      dat <- data.frame(dd)
      # 
      # t <-  try(modbs <- mclustBootstrapLRT(dat$dd, model = "E",
      #                              nboot = 500,
      #                              verbose = FALSE))
      mod <- try(Mclust(dat$dd, G=c(mvbroodmin[i]:mvbroodmax[i]), modelNames = "E"), silent = TRUE)
      if(class(mod) == "try-error") next
      # dat$class <- as.character(mod$classification)
      
      df <- data.frame()
      for (b in 1:length(mod$parameters$mean)){
        df <- rbind(df, data.frame(nsim = ns,
                                   SiteYear = index,
                                   brood = b,
                                   num = length(which(mod$classification == b)),
                                   weight = mod$parameters$pro[b],
                                   mu = mod$parameters$mean[b],
                                   sigma = sqrt(mod$parameters$variance$sigmasq[1]))
        )
      }
      # outlist[[length(outlist) + 1]] <- df
      outdf <- rbind(outdf, df)
    }
  }
  # outdf <- data.table::rbindlist(outlist)
  outdf$filenum <- i
  return(outdf)
}
# }

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
#                           "gdd_tomerge"))
#   
#   test <- parLapply(cl, 1:16, BroodMixMod)
#   stopCluster(cl)
# })
# 
# saveRDS(test, file = "BroodMixSigLim.rds")

test <- readRDS("BroodMixSigLim.rds")
# one question:
# is rpois from GAM prediction equivalent to rbinom from scaled GAM prediction?

a <- rbindlist(test)
rm(test)
a$species <- fs[mvspec[a$filenum]]
a$species <- gsub('gam', '', a$species)
a$species <- gsub('.rds', '', a$species)

# flagging poor clustering models
compound.clusters <- a %>% 
  group_by(species, SiteYear, nsim) %>%
  summarise(mudiff = diff(mu))
zero.cluster <- a %>% 
  filter(num == 0)
# all cases with num = 0 included in compound clusters
flags <- a %>% 
  group_by(species, SiteYear, nsim) %>%
  summarise(mudiff = diff(mu)) %>%
  # filter(mudiff <= 101.9 |
  #        mudiff >= 1142.488)
  filter(mudiff<=20)
a.cut <- anti_join(a, flags, by = c("species", "SiteYear", "nsim"))
#this only seemed to change Horace's DW, maybe bc
#it has really tight brood 3 and 4
#did not clarify cases where ~1/2 the sims split bt extra or not

#try using regionwide G
regions <- site_geo %>%
  select(SiteID, region9, region4)



# b <- a.cut %>% 
b <- a %>%
  group_by(species, SiteYear, nsim) %>%
  mutate(numbrood = max(brood))
b <- b %>%
  group_by(species) %>%
  do( complete(., nesting(SiteYear, nsim, numbrood), brood, 
               fill = list(num = 0, weight = 0, mu = NA, sigma = NA))) %>%
  select(-filenum, -species)

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
  mutate(probsel = nselected / sum(nselected))
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
  summarise(meanN = mean(num),
            meanweight = mean(weight),
            meanmu = mean(mu, na.rm = TRUE),
            meansigma = mean(sigma, na.rm = TRUE),
            sdN = sd(num),
            sdweight = sd(weight),
            sdmu = sd(mu, na.rm = TRUE),
            sdsigma = sd(sigma, na.rm = TRUE),
            prob = probsel[1])

#regional G selection prob
e.reg <- e %>%
  group_by(species, SiteYear) %>%
  filter(numbrood == max(brood)) %>%
  summarise(regprob = prob[1])
t <- str_split_fixed(e.reg$SiteYear, pattern = "_", 2)
e.reg$SiteID <- t[,1]
e.reg$Year <- as.factor(t[,2])
e.reg <- merge(e.reg, regions, by = "SiteID")
e.reg$Reg9Year <- as.factor(paste(e.reg$region9, e.reg$Year, sep = "_"))
e.test <- merge(e.reg, regsites, by = c("Reg9Year"))


f <- e %>%
  group_by(species, SiteYear, brood) %>%
  summarise(N = sum(meanN * prob),
            weight = sum(meanweight * prob, na.rm = TRUE))


t <- str_split_fixed(f$SiteYear, pattern = "_", 2)
f$SiteID <- t[,1]
f$Year <- as.factor(t[,2])

broods <- merge(f, sites, by = "SiteID")

broods <- broods %>%
  group_by(species, SiteYear) %>%
  arrange(brood) %>%
  mutate(lastprop = N[max(brood)] / sum(N[max(brood)-1] + N[max(brood)]),
         lastN = N[max(brood)])

# plot proportion of last brood by species x region
specplot <- ggplot(broods, aes(x = lastprop, group = region4, color = region4)) + geom_density() +
  facet_wrap(~species, scales = "free_y") + theme_bw() +
  xlab("Proportion of last brood size to sum(last 2 broods' size)")
specplot


# plot proportion of last brood 
sitebrood <- broods %>%
  filter(brood == max(brood)) 

broodplot <- ggplot(sitebrood, aes(x = lon, y = lat, color = lastprop)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Year)
broodplot

tt <- e
t <- str_split_fixed(tt$SiteYear, pattern = "_", 2)
tt$SiteID <- t[,1]
tt$Year <- as.factor(t[,2])
rm(t)
tt <- merge(tt, sites, by = "SiteID")
# variable: photoperiod at peaks
library(geosphere)
gdd2vars <- function(SiteYear, meanmu, lat, meanGDD){
  temp <- gdd_tomerge %>% filter(SiteYear == SiteYear)
  ord <- temp$Ordinal[which(temp$cumdegday > meanmu)[1]]
  photo <- daylength(lat, ord)
  gddleft <- max(temp$cumdegday) - meanmu
  gddexpected <- meanGDD - meanmu
  data_frame(ord,
             photo,
             gddleft,
             gddexpected)
}

# tt <- test %>%
#   rowwise() %>%
#   mutate(peakdate = gdd2date(SiteYear, meanmu))

system.time({uu <- tt %>%
  group_by(species, SiteYear, brood, numbrood) %>%
  do(gdd2vars(.$SiteYear, .$meanmu, .$lat, .$meanGDD))
})
saveRDS(uu, "gddvars.rds")

vv <- merge(tt, uu, by = c("species", "SiteYear", "brood", "numbrood"))
# 
# t <- str_split_fixed(tt$SiteYear, pattern = "_", 2)
# tt$SiteID <- t[,1]
# tt$Year <- as.factor(t[,2])
# tt <- merge(tt, sites, by = "SiteID")
# uu <- tt %>%
#   rowwise() %>%
#   mutate(peakphoto = )

broods.cut <- broods %>%
  select(species, SiteYear, brood, N, weight, lastprop, lastN)
vv.cut <- vv %>%
  group_by(species) %>%
  mutate(specmax = max(brood)) %>%
  ungroup() %>%
  rowwise() %>%
  filter(numbrood == specmax)
vv.cut <- vv.cut %>%
  filter(brood == numbrood-1)
moddat <- merge(vv.cut, broods.cut, by = 
                  c("species", "SiteYear", "brood"))

moddat <- moddat %>%
  filter(species == "Wild Indigo Duskywing") %>%
  mutate(zmeanmu = scale(meanmu),
         zphoto = scale(photo),
         zlat = scale(lat))

# model proportion of last brood
library(lme4)
mod <- glmer(cbind(round(lastN), round(N)) ~
               zphoto*zlat +
               (1|region9) +
               (1|Year), 
             data = moddat,
             family = binomial,
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000000))
)
summary(mod)
sjp.int(mod, type = "cond")  
r.squared.merMod(mod)


# model population growth from 1st to 1st broods


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
         lag.siteyear = lag(SiteYear, 1)) %>%
  ungroup() %>%
  mutate(lambda = log((N + 1) / (lag.brood1 + 1))) 

# weird errors here with scale
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# break up by species or model together????

#could scale lastprop by site or region, what makes sense?
newdat <- newdat %>%
  group_by(species, SiteID) %>%
  mutate(sitescaleDD = scale_this(lag.brood1))


#set NA to 0 for DD
newdat$sitescaleDD[is.na(newdat$sitescaleDD)] <- 0
# newdat$scaleProp[is.na(newdat$scaleProp)] <- 0
newdat$scaleprop <- scale_this(newdat$lastprop)
newdat$scalepropsq <- newdat$scaleprop ^ 2
newdat$scalelat <- scale_this(newdat$lat)
vv.gdd <- vv.cut %>%
  select(species, SiteYear, gddleft) %>%
  rename(lag.siteyear = SiteYear)
newdat <- merge(newdat, vv.gdd, by = c("species", "lag.siteyear"))
newdat$scalegddleft <- scale_this(newdat$gddleft)

# newdat <- newdat[complete.cases(newdat),]

popdat <- newdat %>% filter(species == "Northern Broken-Dash")

library(lme4)
mod <- lmer(lambda ~ 
              # sitescaleDD + 
              scaleprop * 
              scalegddleft + 
              (1|Year) + (1|species),
            data = newdat)

mod <- lmer(lambda ~ 
              sitescaleDD +
              scaleprop * 
              scalegddleft + 
              (1|Year),
            data = popdat)

summary(mod)
sjp.int(mod, type = "eff")  
r.squared.merMod(mod)

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
eff <- allEffects(mod)
plot(eff)


library(sjPlot)
# sjp.setTheme(theme = theme_minimal())

sjt.lmer(mod1, file = "brood.html")


sjp.lmer(mod1,facet.grid = FALSE,
         sort.coef = "sort.all")
# sjp.lmer(weather, "coef")

sjp.lmer(mod1, type = "poly", poly.term = "scaleprop")
sjp.int(mod, type = "cond")  


# not much added with gam
gmod <- gamm(lambda ~ s(scaleDD) + 
               s(scaleProp) + 
               s(lat), random = list(Year=~1),
             data = newdat)


