source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
library(pastecs)
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
for (i in 1:nrow(species)){
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
  
  # trying to add GDD instead of ordinal date
  counts <- merge(counts, gdd, by = c("SiteID", "SiteDate"), all.x = TRUE, all.y = FALSE)
  
  print(sp)
  
  counts <- counts[, `:=` (SurvPerYear = length(unique(SeqID)),
                           YearTotal = sum(Total)), 
                   by = list(SiteID, Year)]
  
  datGAM <- counts[YearTotal >= 1]
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
      filter(YearTotal > 0,
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
    tempdf$cumdegday <- c(rep(0, nrow(tempdf)/2), rep(max(temp$cumdegday), nrow(tempdf)/2))
    
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
  
  
  temp <- dat
  
  if(sum(temp$Total) < 20) next
  if(length(unique(temp$SiteID)) < 2) next
  if(length(unique(temp$Year)) < 2) next
  
  mod7 <- try(gam(Total ~ #s(listlength, k = 3) + #s(temperature, k = 3) + s(duration, k = 3) +
                    # s(SiteID, bs = "re") + s(Year, bs = "re") +
                    s(SiteYear, Ordinal, bs = "fs", k = 10, m = 1) +
                    te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(7, 20), d = c(2, 1)),
                  family = nb(theta = NULL, link = "log"),
                  # family = poisson(link = "log"),
                  data = temp,
                  method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
  # 
  # model used for enormous gamlist results
  # gam(Total ~ s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
  #   te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
  # family = nb(theta = NULL, link = "log"),
  # # family = poisson(link = "log"),
  # data = temp,
  # method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
  # 
  # 
  # # GDD instead of ordinal is slightly worse for Giant Swallowtail
  # mod8 <- try(gam(Total ~ s(listlength, k = 3) +  #s(temperature, k = 3) + s(duration, k = 3) +
  #                   s(SiteID, bs = "re") + s(Year, bs = "re") +
  #                   s(Year, cumdegday, bs = "fs", m = 1) +
  #                   te(lat, lon, cumdegday, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
  #                 family = nb(theta = NULL, link = "log"),
  #                 # family = poisson(link = "log"),
  #                 data = temp,
  #                 method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
  # 
  # 
  if(("gam" %in% class(mod7)) == FALSE){
    mod7 <- try(mod7 <- gam(Total ~ #s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                              s(SiteYear, Ordinal, bs = "fs", k = 5, m = 1) +
                              te(lat, lon, Ordinal, bs = c("tp", "tp"), k = c(5, 15), d = c(2, 1)),
                            family = nb(theta = NULL, link = "log"),
                    # family = poisson(link = "log"),
                    data = temp,
                    method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
  }

  if(("gam" %in% class(mod7)) == FALSE){
    mod7 <- try(mod7 <- gam(Total ~ #s(listlength, k = 3) + s(temperature, k = 3) + s(duration, k = 3) +
                              SiteID + Year +
                      s(Ordinal, bs = "tp", k = 15),
                    family = nb(theta = NULL, link = "log"),
                    # family = poisson(link = "log"),
                    data = temp,
                    method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500)))
  }

  if(("gam" %in% class(mod7)) == FALSE) next
  
  temp <- temp %>% filter(Ordinal > 70 & Ordinal < 320)
  
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
  
  newData$GAM.pred <- as.vector(predict.gam(mod7, newData, type = "response"))
  # temp$GAM.pred <- predict.gam(mod7, type = "response")
  # newData$GAM.pred <- predict.gam(mod7, newData, type = "response", 
  #                                 exclude = c("te(lat,lon,Ordinal"))
  newData <- data.table(newData)
  newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteYear"] 
  
  c <- ggplot(data = newData, aes(x = Ordinal, y = Gamma, group = SiteID, color = SiteID)) +
    geom_line(size = 1, alpha = .5) + 
    theme_bw() + theme(legend.position = "none") +
    facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", sp, sep = " "))
  # print(c)
  d <- ggplot(data = newData, aes(x = Ordinal, y = GAM.pred, group = SiteID, color = SiteID)) + 
    geom_line(size = 1, alpha = .5) +
    theme_bw() + geom_point(data = temp, aes(x = Ordinal, y = Total))
  d <- d + facet_wrap( ~ Year, scales = "free_y") + theme(legend.position = "none") +
    ggtitle(paste("GAM Predictions", sp, sep = " "))
  # print(d)
  
  pdf(paste("CovGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(c)
  dev.off()
  
  pdf(paste("CovGAMPRED", sp, ".pdf", sep = ""), width = 10, height = 6)
  print(d)
  dev.off()
  
  # turning points
  
  t <- lapply(unique(newData$SiteYear), tpGAM, Dat = newData)
  tt <- data.table(rbindlist(t))
  
  gamlist$mod <- mod7
  gamlist$preds <- newData
  gamlist$counts <- temp
  gamlist$turning <- tt
  
  saveRDS(gamlist, paste("gam", sp, ".rds", sep = ""))
}

# saveRDS(gamlist, "gamlist_withdetcov.rds")

# for gamlist, half size of file by removing redundant gam model

fs <- list.files('gamGDD')
for (i in 1:length(fs)){
  gamlist <- readRDS(paste("gamGDD", fs[i], sep = "/"))
  if("gam" %in% class(gamlist$mod)){
    if("gam" %in% class(gamlist$modgdd)){
      gamlist$modgdd <- NULL
      saveRDS(gamlist, file = paste("gamGDD", fs[i], sep = "/"))
    }
  }
}




gamlist <- readRDS("gamlist_nodetcov.rds")
for (i in 1:length(gamlist)){
  templist <- gamlist[[i]]
  if(length(templist) > 1){
  sp <- templist[[1]]
  # if("gam" %in% class(templist$mod)){
  #   if("gam" %in% class(templist$modgdd)){
  #     gamlist$modgdd <- NULL
      saveRDS(templist, file = paste("gamGEOnodet/", sp, ".rds", sep = ""))
  #   }
  # }
  }
}


# plot number of peaks, timing of peaks, weights with and without cutoff threshold for weight
# what threshold for a real peak works?
# i <- 96
# temp <- gamlist[[i]]$turning
temp <- tt
print(species$CommonName[i])

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



# test with Giant Swallowtail population model
temp <- tt
t <- str_split_fixed(temp$SiteYear, pattern = "_", 2)
temp$SiteID <- t[,1]
temp$Year <- t[,2]

temp$SiteDate <- as.Date(strptime(paste(temp$Year, temp$Ordinal, sep = " "), format = "%Y %j"))

# remaining gdd left in season?
gdd$date <- NULL
gdd <- gdd %>%
  group_by(SiteID, year) %>%
  mutate(remaindegday = max(cumdegday) - cumdegday) %>%
  select(SiteDate, SiteID, cumdegday, remaindegday)

temp <- merge(temp, gdd, by = c("SiteID", "SiteDate"), all.x = TRUE, all.y = FALSE)
temp$year <- NULL

temp <- temp %>%
  filter(tptype == "peak") %>%
  group_by(SiteID, Year) %>%
  mutate(broodnum = 1:length(tptype)) %>%
  select(SiteID, Year, SiteYear, aoc, weight, cumdegday, remaindegday, broodnum)
  
maxBrood <- max(temp$broodnum)
results <- list()
for (i in 1:length(unique(temp$SiteYear))){
  sy <- unique(temp$SiteYear)[i]
  t <- temp %>% filter(SiteYear == sy)
  if(nrow(t) < maxBrood){
    newrow <- t[1,]
    newrow$aoc <- 0
    newrow$weight <- 0
    newrow$cumdegday <- NA
    newrow$remaindegday <- NA
    newrow$broodnum <- nrow(t) + 1
    t <- rbind(t, newrow)
  }
  # t$broodnum <- paste("brood", t$broodnum, sep = "")
  nextbrood <- temp %>% filter(SiteID == t$SiteID[1],
                               Year == as.character(as.numeric(t$Year[1]) + 1))
  if(nrow(nextbrood) == 0) next
  
  lmdat <- nextbrood %>% filter(broodnum == 1)
  for (j in 1:maxBrood){
    tcol <- t %>% filter(broodnum == j) %>%
      ungroup() %>%
      select(aoc, weight, cumdegday, remaindegday)
    names(tcol) <- paste(names(tcol), j, sep = "_")
    lmdat <- cbind(lmdat, tcol)
  }
  results[[length(results) + 1]] <- lmdat
}
moddat <- rbindlist(results)
moddat$cumdegday_3[which(is.na(moddat$cumdegday_3) == TRUE)] <- mean(moddat$cumdegday_3, na.rm = TRUE)
moddat$remaindegday_3[which(is.na(moddat$remaindegday_3) == TRUE)] <- mean(moddat$remaindegday_3, na.rm = TRUE)
moddat$lambda <- log((moddat$aoc + 1)/(moddat$aoc_1 + 1))

library(lme4)
mod <- lmer(lambda ~ log(aoc_1+1) + remaindegday_2  * log(aoc_2+1) + 
            remaindegday_3 * log(aoc_3+1) +
            cumdegday + (1|Year), data = moddat)




# get some GAM predictions
gamlist <- readRDS("gamGDD/gamRed-spotted Purple.rds")
# gamlist <- readRDS("gamGEOnodet/Red-spotted Purple.rds")
pred <- as.data.frame(gamlist$preds)
datGAM <- as.data.frame(gamlist$counts)
####
newData <- pred
c <- ggplot(data = newData, aes(x = Ordinal, y = Gamma, group = SiteID, color = SiteID)) +
  geom_line(size = 1, alpha = .5) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", gamlist[[1]], sep = " "))
print(c)
c <- ggplot(data = newData, aes(x = cumdegday, y = Gamma, group = SiteID, color = SiteID)) +
  geom_line(size = 1, alpha = .5) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", gamlist[[1]], sep = " "))
print(c)

#####################################################
# test how GAM predictions do when thrown into stopover model
counts <- pred %>%
  filter(Year == "2003") %>%
  data.frame()
# counts <- datGAM %>%
#   filter(Year == "2003") %>%
#   data.frame()
# see how survey intervals relate to gdd intervals
# t <- counts %>% 
#   group_by(SiteID) %>% 
#   arrange(Ordinal) %>% 
#   summarise(meandiffday = mean(diff(Ordinal)), 
#             meandiffgdd = mean(diff(cumdegday)))


counts$SiteDate <- strptime(paste(counts$Year, counts$Ordinal, sep = " "), format = "%Y %j")
counts$SiteDate <- as.Date(counts$SiteDate)
counts$Week <- week(counts$SiteDate)
counts$Week <- counts$Week - min(counts$Week) + 1
#overachieving volunteers going out more than once a week!
#choose first one, not averaging, easier to match with covariates
counts_uniq <- counts %>% 
  group_by(SiteID, Week) %>%
  arrange(SiteDate) %>%
  summarise(Total = rpois(1, GAM.pred[1]))
  # summarise(Total = Total[1])

count_matrix <- as.matrix(cast(counts_uniq, SiteID ~ Week, value = "Total"))
# count_matrix <- round(count_matrix)
# count_matrix[sample(1:length(count_matrix), .2*length(count_matrix))] <- NA

# count_matrix[is.na(count_matrix)] <- NA

site_covs <- counts %>%
  select(SiteID, lat) %>%
  distinct()

counts <- count_matrix

cov.p <- array(rnorm(dim(counts)[1]*dim(counts)[2]), dim = c(dim(counts), 1))
cov.p[,,1][is.na(counts)] <- NA
# select sites with enough individuals counted
siteRows <- which(rowSums(counts, na.rm = TRUE) >= 5)
counts <- counts[siteRows, ]
counts[is.na(counts)] <- -1


M <- 3
S <- dim(counts)[1]
K <- TIME <- dim(counts)[2]

# Covariate weather for p
# cov.p <- matrix(data_samp$TEMP,nrow=S,ncol=K,byrow=TRUE)
#cov.p now needs to be a S by T by qp (number of covariates) array
# covs <- c(2,4) #Select detection covariate here (1:5 possible)
# covs <- c(1,7) #temperature and list-length of species observed
# if (length(covs) > 1) cov.p <- cov_array[,,covs]
# if (length(covs) == 1) cov.p <- array(data = cov_array[,,covs], dim = c(dim(cov_array[,,covs]), 1))
cov.p <- cov.p[siteRows, , , drop = FALSE]
# qp is the number of covariates for p
qp <- dim(cov.p)[3]
for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
cov.p[is.na(cov.p)] <- -1


# Time covariate for phi
cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))),S,K-1,byrow=TRUE) 

# Covariate (latitude)
# lat <- rowMeans(cov_array[,,5], na.rm = TRUE)
lat <- site_covs[, "lat"]
lat <- lat[siteRows]
cov.w <- cov.mu <- scale(lat)[,1]


# add bootstrap
# 
# transCounts = counts
# trans.cov.phi = cov.phi
# trans.cov.w = cov.w
# trans.cov.mu = cov.mu
# trans.cov.p = cov.p
# 
# 
# nBoots = 1
# bootFits = list()
# bootIdxs = array(0,c(nBoots,dim(counts)[1]))
# bootTime <- list()
# # 
# # bootFits <- readRDS("GSFBoot.rds")
# # bootIdxs <- readRDS("GSFIDxs.rds")
# # bootTime <- readRDS("GSFTime.rds")
# 
# 
# # 8 bootstrap replicates took 2 days!
# for (b in 48:nBoots) {
#   
#   print(b)
#   startTime <- Sys.time()
#   
#   bSites = sample(dim(counts)[1],replace=T)
#   bootIdxs[b,] = bSites
#   counts = transCounts[bSites,]
#   cov.w = trans.cov.w[bSites]
#   cov.mu = trans.cov.mu[bSites]
#   cov.p = trans.cov.p[bSites,,, drop = FALSE]
#   cov.phi = trans.cov.phi[bSites, ]

## MODEL RUNNING ######################################
# source('../../Dropbox/SESYNC ZGroup Summer2015/Ohio/FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

########################################################
# p.m <- "cov"
p.m <- "common"
w.m <- "cov"
mu.m <- "cov"
sigma.m <- "het"
# phi.m <- "logit.a"
phi.m <- "const"

Tries <- 2
temp.fit <- list()
temp.ll <- rep(NA, Tries)

for (k in 1:Tries){
  start.list <- startVals(p.m,w.m,mu.m,sigma.m,phi.m)
  pars.start <- c(start.list$N,start.list$cvec, start.list$d0, start.list$d1, start.list$b0, start.list$b1,start.list$sigma,  start.list$a0, start.list$a1, start.list$a2) #this line remains the same for all models
  temp.fit[[k]] <- mLLMixtCounts.fit(p.m,w.m,mu.m,sigma.m,phi.m)
  temp.ll[k] <- temp.fit[[k]]$ll.val
}

if (length(which(is.na(temp.ll) == TRUE)) < Tries){
  tempchoose <- min(c(1:Tries)[which(temp.ll==max(temp.ll, na.rm=TRUE))])
  temp <-temp.fit[[tempchoose]]
}else{
  temp <- list()
  temp$ll.val <- NA
}

tempchoose <- min(c(1:5)[which(temp.ll==max(temp.ll,na.rm=T))])

temp <-temp.fit[[tempchoose]] 
################################################################
###############################################################
library(mclust)
library(plyr)
library(dplyr)
library(mgcv)
library(lubridate)
library(tidyr)
library(stringr)
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
# a <- ggplot(data = sites, aes(x = lon, y = lat, group = class, color = class)) + geom_point()
# a

library(data.table)
site_geo <- read.csv("C:/Users/Tyson/Desktop/Box Sync/Ohio/GIS/OHsites_reconciled.csv", header = TRUE)
site_geo <- data.table(site_geo)
setnames(site_geo,"Name","SiteID")
site_geo[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

site_geo <- merge(site_geo, sites, by = c("SiteID", "Description.x", "lat", "lon"))




####
# clustering broods from GAM smooths
####

# 2 things
# might could separate giant gamlist into species and resave so more accessible
# saved model twice in each species gam, lots more disk space than necessary
preddat <- pred %>% filter(Year == 2003)
a <- ggplot(data = preddat, aes(x = cumdegday, y = GAM.pred.reg9yr, group = SiteID, color = SiteID))
a + geom_line() + facet_wrap(~ Reg9Year)


# ScaleDennis <- function(x){exp(x)/sum(exp(x))}
# pred <- pred[, Gamma2 := ScaleDennis(GAM.pred), by = "SiteYear"] 


# another issue, gdd alters scale of x-axis in weird ways
# use "butterfly-days" or "gdd-days" as population index?
# can use approximate merge to get gdd and ordinal date to match up
# how to use this in mixture model interpretation?

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
                                      # exclude = c("s(SiteID)"))
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

#RSP try, mclust with limited 2:3 G, sigma same, gdd
system.time({dun1 <- BroodMixMod(2)})




library(parallel)


# multicore
system.time({
cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(mclust)
  library(data.table)
  library(dplyr)
  library(mgcv)
})
clusterExport(cl=cl, 
              varlist=c("mvspec", 
                        "mvbroodmax", 
                        "mvbroodmin", 
                        "gdd_tomerge"))

test <- parLapply(cl, 1:16, BroodMixMod)
stopCluster(cl)
})

saveRDS(test, file = "BroodMixSigLim.rds")

test <- readRDS("BroodMixSigLim.rds")
# one question:
# is rpois from GAM prediction equivalent to rbinom from scaled GAM prediction?

# compare RSP broods with and without mclust limited to 2 or 3 broods
saveRDS(outdf, file = "RSPlimited.rds")
saveRDS(outdf, file = "RSP1to5.rds")
saveRDS(outdf, file = "RSPlimited_sigmahet.rds")



outdf <- test[[16]]
# best way to summarise mixture model simulations might be conditional probability
# mean weight of a brood, conditional on it being selected as a mode in mclust

a <- outdf %>% filter(SiteYear == "001_1997") %>% data.frame()

b <- a %>% 
  group_by(nsim) %>%
  mutate(numbrood = max(brood)) %>%
  ungroup() %>%
  complete(nsim, brood, fill = list(num = 0, weight = 0, mu = NA, sigma = NA))


c <- b %>%
  group_by(nsim) %>%
  summarise(maxbrood = numbrood[1]) %>%
  group_by(maxbrood) %>%
  summarise(nselected = n())
c$probsel <- c$nselected / sum(c$nselected)
c <- dplyr::rename(c, numbrood = maxbrood)

d <- merge(b, c, by = "numbrood")
  
e <- d %>%
  group_by(brood, numbrood) %>%
  summarise(meanN = mean(num),
            meanweight = mean(weight),
            meanmu = mean(mu, na.rm = TRUE),
            meansigma = mean(sigma, na.rm = TRUE),
            medianN = median(num),
            medianweight = median(weight),
            medianmu = median(mu, na.rm = TRUE),
            mediansigma = median(sigma, na.rm = TRUE),
            prob = probsel[1])

f <- e %>%
  group_by(brood) %>%
  summarise(N = sum(meanN * prob),
            weight = sum(meanweight * prob, na.rm = TRUE))



# compare different mclust estimates for rsp

# compare RSP broods with and without mclust limited to 2 or 3 broods
outdf1 <- readRDS(file = "RSPlimited.rds")
outdf1$modelrun <- "sigmahom_2to3"
outdf2 <- readRDS(file = "RSP1to5.rds")
outdf2$modelrun <- "sigmahom_1to5"
outdf3 <- readRDS(file = "RSPlimited_sigmahet.rds")
outdf3$modelrun <- "sigmahet_2to3"

outdf <- rbind(outdf1, outdf2, outdf3)

a <- outdf %>% 
  group_by(modelrun, SiteYear, nsim) %>%
  mutate(numbrood = max(brood)) 
  # ungroup() %>%
  # complete(nsim, brood, fill = list(num = 0, weight = 0, mu = NA, sigma = NA))

b <- a %>%
  group_by(modelrun, SiteYear, brood, numbrood) %>%
  summarise(meanN = mean(num),
            meanweight = mean(weight),
            meanmu = mean(mu, na.rm = TRUE),
            meansigma = mean(sigma, na.rm = TRUE),
            nselect = n())
c <- b %>%
  group_by(modelrun, SiteYear, brood) %>%
  summarise(N = sum(meanN * nselect / 100),
            weight = sum(meanweight * nselect / 100, na.rm = TRUE))

c <- b %>%
  group_by(modelrun, brood) %>%
  summarise(
            weight = mean(meanweight * nselect / 100, na.rm = TRUE))

# hom sigma, limited brood number seems best


saveRDS(outdf, "leastskip.limit.rds")
saveRDS(outdf, "peckskip.limit.rds")
# plot composition of generations across space and time
outdf1 <- readRDS("leastskip.limit.rds")


outdf1 <- test[[1]]
a <- outdf1

a <- rbindlist(test)
a$species <- fs[mvspec[a$filenum]]
a$species <- gsub('gam', '', a$species)
a$species <- gsub('.rds', '', a$species)


b <- a %>% 
  group_by(species, SiteYear, nsim) %>%
  mutate(numbrood = max(brood)) %>%
  ungroup() %>%
  group_by(species) %>%
  complete(species, SiteYear, nsim, brood, fill = list(num = 0, weight = 0, mu = NA, sigma = NA, numbrood = max(numbrood)))
# problem here, not completing, leaving euro skip with 1 brood, messes up calcs below

c <- b %>%
  group_by(species, SiteYear, nsim) %>%
  summarise(maxbrood = numbrood[1]) %>%
  group_by(species, SiteYear, maxbrood) %>%
  summarise(nselected = n()) %>%
  group_by(species, SiteYear) %>%
  mutate(probsel = nselected / sum(nselected))
c <- dplyr::rename(c, numbrood = maxbrood)

d <- merge(b, c, by = c("species", "SiteYear", "numbrood"))

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
  mutate(lastprop = N[max(brood)] / sum(N[max(brood)-1] + N[max(brood)]))

specplot <- ggplot(broods, aes(x = lastprop, group = region4, color = region4)) + geom_density() +
  facet_wrap(~species, scales = "free_y")
specplot


# plot proportion of last brood
sitebrood <- broods %>%
  filter(brood == max(brood)) 

broodplot <- ggplot(sitebrood, aes(x = lon, y = lat, color = lastprop)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Year)
broodplot


# model population growth from 1st to 1st broods


# expanddat <- broods %>%
#   tidyr::expand(SiteID, Year) %>%
#   dplyr::left_join(broods)

newdat <- 
  broods %>%
  group_by(SiteYear) %>%
  mutate(siteyearN = sum(N)) %>%
  ungroup %>%
  filter(brood == 1) %>%
  arrange(SiteID, Year) %>%
  group_by(SiteID) %>%
  mutate(lag.brood1 = lag(N, 1)) %>%
  ungroup() %>%
  mutate(lambda = log((N + 1) / (lag.brood1 + 1))) 

# weird errors here with scale
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

#could scale lastprop by site or region, what makes sense?
newdat <- newdat %>%
  group_by(SiteID) %>%
  mutate(sitescaleDD = scale_this(lag.brood1))
#set NA to 0 for DD
newdat$sitescaleDD[is.na(newdat$sitescaleDD)] <- 0
# newdat$scaleProp[is.na(newdat$scaleProp)] <- 0
newdat$scaleprop <- scale_this(newdat$lastprop)
newdat$scalepropsq <- newdat$scaleprop ^ 2
newdat$scalelat <- scale_this(newdat$lat)
newdat <- newdat[complete.cases(newdat),]

library(lme4)
mod <- lmer(lambda ~ sitescaleDD + 
              scaleprop * 
              scale(lat) + 
              (1|Year),
            data = newdat)

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
eff <- allEffects(mod1)
plot(eff)


library(sjPlot)
# sjp.setTheme(theme = theme_minimal())

sjt.lmer(mod1, file = "brood.html")


sjp.lmer(mod1,facet.grid = FALSE,
         sort.coef = "sort.all")
# sjp.lmer(weather, "coef")

sjp.lmer(mod1, type = "poly", poly.term = "scaleprop")
sjp.int(mod1, type = "eff")  


# not much added with gam
gmod <- gamm(lambda ~ s(scaleDD) + 
              s(scaleProp) + 
              s(lat), random = list(Year=~1),
            data = newdat)




# Fake data
k=1000
set.seed(104)
val = rnorm(k)
dens = density(val, n=512)

# Convert to counts
dens$y = k/sum(dens$y) * dens$y

plot(dens)


library(mixsmsn)


mixer <- smsn.mix(ddcount, nu = 10, g = 2, get.init = TRUE,
      criteria = TRUE, group = TRUE, family = "Skew.normal", 
      calc.im = FALSE)
mix.hist(ddcount, mixer)
mixauto <- smsn.search(ddcount, nu = -5, g.min = 1, g.max = 5, 
            family = "Skew.normal", criteria = "bic")
mix.hist(ddcount, mixauto$best.model)

library(distr)

## Construct the distribution object.
myMix <- UnivarMixingDistribution(Norm(mean=642, sd=115), 
                                  Norm(mean=977, sd=111),
                                  Norm(mean=1488, sd=265),
                                  mixCoeff=c(0.3, 0.56, 0.14))
myMix <- UnivarMixingDistribution(Norm(mean=888, sd=212),
                                  Norm(mean=1610.5, sd=181),
                                  mixCoeff = c(.91, .09))
## ... and then a function for sampling random variates from it
rmyMix <- r(myMix)

## Sample a million random variates, and plot (part of) their histogram
x <- rmyMix(1e6)
hist(x, breaks=100, col="grey", main="")






library(mixtools)
test <- boot.comp(ddcount, B = 100, max.comp=5,mix.type="normalmix",
                  maxit=400,epsilon=1e-2, arbvar = FALSE)
test2 <- boot.comp(ddcount, B = 100, max.comp=5,mix.type="normalmix",
                  maxit=400,epsilon=1e-2, arbvar = TRUE)
test2 <- normalmixEM(ddcount, arbvar = FALSE, k = 3)
test3 <- normalmixEM(ddcount, k = 3)


     t <-  try(modbs <- mclustBootstrapLRT(ddcount, model = "E",
                                  nboot = 100,
                                  verbose = FALSE))

test4 <- Mclust(daycount, G = c(1:4), modelNames = "E")

library(flexmix)
pr <- pred %>%
  filter(Year == 2003) %>%
  select(SiteID, Reg9Year, Gamma, GAM.pred.reg9yr, cumdegday) %>% 
  group_by(SiteID) %>%
  mutate(GamTot = round(sum(GAM.pred.reg9yr))) %>%
  rowwise() %>%
  mutate(SimCount = rbinom(1, GamTot, Gamma))

pr <- pred %>%
  filter(Year == 2003) %>%
  select(Reg9Year, Gamma, GAM.pred.reg9yr, cumdegday) %>% 
  distinct() %>%
  rowwise() %>%
  mutate(SimCount = rpois(1, GAM.pred.reg9yr))

region <- rep(pr$Reg9Year, pr$SimCount)
dd <- rep(pr$cumdegday, pr$SimCount)
df <- data.frame(dd, region)

test <- stepFlexmix(dd ~ 1|region,
                model = FLXMRglm(family = "gaussian"), k = 3, nrep = 3,
                data = df)


tapply(datGAM$Reg9Year, datGAM$Reg9Year, length)

index <- "Columbus_2006"

# #clustering with raw counts leads to too many modes
countdat <- datGAM %>%
  filter(Reg9Year == index) %>%
  select(Total, SiteID, Ordinal, lat, lon, cumdegday)
# yd <- rep(clustdat$Ordinal, clustdat$Total)
# lat <- rep(clustdat$lat, clustdat$Total)
# lon <- rep(clustdat$lon, clustdat$Total)
# dd <- rep(clustdat$cumdegday, clustdat$Total)
# 
# dat <- data.frame(cbind(yd, lat, lon, dd))
# mod <- Mclust(dat$dd, G = 1:5)
# summary(mod, parameters = TRUE)
# dat$class <- as.factor(mod$classification)
# dat$uncert <- mod$uncertainty
# d <- ggplot(aes())

# cluster with GAM predictions of region
clustdat <- pred %>%
  filter(Reg9Year == index) %>%
  # select(Ordinal, GAM.pred, lat, SiteID) %>%
  select(cumdegday, GAM.pred.reg9yr) %>%
  distinct()
  # rowwise() %>%
  # mutate(GAMpois = rpois(1, GAM.pred)) %>%
  # ungroup()

clustdat$SiteID <- as.factor(clustdat$SiteID)
a <- ggplot(data = clustdat, aes(x = cumdegday, y = GAM.pred.nosite, group = SiteID, color = SiteID))
a + geom_line() + geom_point(data = countdat, aes(x = cumdegday, y = Total, group = SiteID, color = SiteID))


# gddfilt <- gdd %>%
#   filter(year == 2003) %>%
#   select(SiteID, yday, cumdegday)
# gddfilt$Ordinal <- gddfilt$yday
# clustdat <- merge(clustdat, gddfilt, by = c("SiteID", "Ordinal"))
# yd <- rep(clustdat$Ordinal, clustdat$GAMpois)
# clustdat <- clustdat %>% filter(SiteID == "017")
# lat <- rep(clustdat$lat, clustdat$GAMpois)
# dd <- rep(clustdat$cumdegday, clustdat$GAMpois)
# lat <- rep(clustdat$lat, clustdat$GAM.pred.reg9yr * 10)
dd <- rep(clustdat$cumdegday, clustdat$GAM.pred.reg9yr * 10)

# dat <- data.frame(cbind(yd, lat, dd))
# dat <- data.frame(cbind(lat, dd))
dat <- data.frame(dd)

mod <- Mclust(dd, G = 1:5)
mod <- densityMclust(dat[,2], G = 1:5, modelNames = "V")
dat$class <- as.character(mod$classification)

a <- ggplot(dat, aes(x = lat, y = dd, color = class)) + geom_jitter()
a

b <- ggplot(dat, aes(x = dd, color = class)) + geom_freqpoly()
b <- ggplot(dat, aes(x = dd, fill = class)) + geom_histogram(binwidth = 10)
b + geom_line(data = clustdat, aes(x = cumdegday, y = GAM.pred, group = SiteID, color = SiteID))
b



clustdat <- datGAM %>%
  filter(Total > 0, SiteID == "002", Year == 2001) %>%
  select(Total, yday, lat, lon, cumdegday)
yd <- rep(clustdat$yday, clustdat$Total)
lat <- rep(clustdat$lat, clustdat$Total)
dd <- rep(clustdat$cumdegday, clustdat$Total)


dat <- data.frame(cbind(yd, lat, dd))

mod <- densityMclust(dat[,3], G = 1:4, modelNames = "E")
dat$class <- as.character(mod$classification)

a <- ggplot(dat, aes(x = yd, y = dd, color = class)) + geom_point()
a


library(distr)

## Construct the distribution object.
myMix <- UnivarMixingDistribution(Norm(mean=355, sd=sqrt(8288)), 
                                  Norm(mean=997, sd=sqrt(21061)),
                                  Norm(mean=1508, sd=sqrt(24753)),
                                  mixCoeff=c(0.08, 0.6, 0.32))
## ... and then a function for sampling random variates from it
rmyMix <- r(myMix)

## Sample a million random variates, and plot (part of) their histogram
x <- rmyMix(1e6)
hist(x, breaks=100, col="grey", main="")

