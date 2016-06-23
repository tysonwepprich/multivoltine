source('bootstrapMfunctions.R')

library(mgcv)
library(ggplot2)
library(lubridate)
library(stringr)
library(pastecs)
ScaleSumTo1 <- function(x){x/sum(x)}

fs <- list.files('gamGDD')
species <- str_split_fixed(string = fs, pattern = "gam", 2)[,2]
species <- str_split_fixed(string = species, pattern = ".rds$", 2)[,1]


# plot all species Gamma's to assess number of generations
for (i in 1:length(species)){
sp <- species[i]
# get some GAM predictions
gamlist <- readRDS(paste("gamGDD/gam", sp, ".rds", sep = ""))
# gamlist <- readRDS("gamGEOnodet/Red-spotted Purple.rds")
pred <- as.data.frame(gamlist$preds)
datGAM <- as.data.frame(gamlist$counts)
# filtcount <- datGAM %>%
#   group_by(SiteYear) %>%
#   mutate(SiteYearTotal = sum(Total),
#          SiteYearSurv = length(Total))
####
#plot by degree day or date
newData <- pred
# c <- ggplot(data = newData, aes(x = Ordinal, y = Gamma, group = SiteID, color = SiteID)) +
#   geom_line(size = 1, alpha = .5) + 
#   theme_bw() + theme(legend.position = "none") +
#   facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", gamlist[[1]], sep = " "))
# print(c)
c <- ggplot(data = newData, aes(x = cumdegday, y = Gamma, group = SiteID, color = lat)) +
  geom_line(size = 1, alpha = .5) + 
  theme_bw() + theme(legend.position = "none") +
  facet_wrap( ~ Year, scales = "free_y") + ggtitle(paste("Scaled Phenology", gamlist[[1]], sep = " "))
# print(c)


pdf(paste("CovGamma", sp, ".pdf", sep = ""), width = 10, height = 6)
print(c)
dev.off()
}


# 1st pass of species, send to sesync cluster 5/9/16
specs <- species[sort(c(22, 38, 74, 69, 59, 52, 49, 43))]

firstrun <- as.list(specs)
for (i in 1:length(specs)){
  sp <- firstrun[[i]]
  # get some GAM predictions
  gamlist <- readRDS(paste("gamGDD/gam", sp, ".rds", sep = ""))
  # gamlist <- readRDS("gamGEOnodet/Red-spotted Purple.rds")
  pred <- as.data.frame(gamlist$preds)
  firstrun[[i]]$GAMpred <- pred
}
saveRDS(firstrun, "firstrundata.rds")

specBrood <- data.frame(species = specs, 
                        minBrood = c(1, 2, 1, 1, 1, 2, 2, 3),
                        maxBrood = c(3, 4, 2, 2, 3, 3, 3, 4))
saveRDS(specBrood, "firstrunbrood.rds")





# 2nd pass of species, send to sesync cluster 5/9/16
specs <- species[sort(c(10, 27, 46, 49, 62, 75, 81))]

secondrun <- as.list(specs)
for (i in 1:length(specs)){
  sp <- secondrun[[i]]
  # get some GAM predictions
  gamlist <- readRDS(paste("gamGDD/gam", sp, ".rds", sep = ""))
  # gamlist <- readRDS("gamGEOnodet/Red-spotted Purple.rds")
  pred <- as.data.frame(gamlist$preds)
  secondrun[[i]]$GAMpred <- pred
}
saveRDS(secondrun, "secondrundata.rds")

specBrood <- data.frame(species = specs, 
                        minBrood = c(2, 2, 2, 1, 2, 2, 2),
                        maxBrood = c(3, 3, 3, 2, 3, 3, 3))
saveRDS(specBrood, "secondrunbrood.rds")










# from GAM plots, estimate M for each year to constrain stopover model
params <- data.frame(yr = as.character(1996:2014), 
                     minBroods = rep(3, 19), 
                     nsim = 5)

# data_file Rdata
dataIN <- c("pred", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

# # single core
system.time({
  test <- lapply(paramIN$nRun, StopoverGAM)
  })


# multicore
system.time({
cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(reshape)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, StopoverGAM)
stopCluster(cl)
})


# function to fit stopover model to GAM predictions
# Year is character, obsM best guess from GAM plots
StopoverGAM <- function(nRun){
  pars <- params[nRun, ]
  nsim <- pars$nsim
  obsM <- pars$obsM
  yr <- pars$yr
  outlist <- as.list(1:nsim)
  for (i in 1:nsim){
    counts <- pred %>%
      filter(Year == yr) %>%
      data.frame()
    
    # if modeled by gdd
    counts <- counts[which(counts$cumdegday %% 50 == 0), ]
    counts$Week <- counts$cumdegday / 50 + 1
    
    counts_uniq <- counts %>% 
      group_by(SiteID, Week) %>%
      arrange(Week) %>%
      summarise(Total = rpois(1, GAM.pred[1])) # here's the random part in each sim
    
    count_matrix <- as.matrix(reshape::cast(counts_uniq, SiteID ~ Week, value = "Total"))
    
    site_covs <- counts %>%
      select(SiteID, lat) %>%
      distinct()
    
    counts <- count_matrix
    
    # dummy p covariate to make model work, but really just want common Np index
    cov.p <- array(rnorm(dim(counts)[1]*dim(counts)[2]), dim = c(dim(counts), 1))
    cov.p[,,1][is.na(counts)] <- NA
    # select sites with enough individuals counted
    siteRows <- which(rowSums(counts, na.rm = TRUE) >= 5)
    counts <- counts[siteRows, ]
    counts[is.na(counts)] <- -1
    
    M <<- obsM
    S <<- dim(counts)[1]
    K <<- dim(counts)[2]
    TIME <<- dim(counts)[2]
    
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
    
    # don't know why double assignment makes things work in parallel
    counts <<- counts
    cov.w <<- cov.w
    cov.mu <<- cov.mu
    cov.p <<- cov.p
    cov.phi <<- cov.phi
    qp <<- dim(cov.p)[3]
    
    
    # p.m <- "cov"
    p.m <<- "common"
    w.m <<- "cov"
    mu.m <<- "cov"
    sigma.m <<- "het"
    # phi.m <- "logit.a"
    phi.m <<- "const"
    
    Tries <- 3
    temp.fit <- list()
    temp.ll <- rep(NA, Tries)
    
    for (k in 1:Tries){
      start.list <<- startVals(p.m,w.m,mu.m,sigma.m,phi.m)
      pars.start <<- c(start.list$N,start.list$cvec, start.list$d0, start.list$d1, start.list$b0, start.list$b1,start.list$sigma,  start.list$a0, start.list$a1, start.list$a2) #this line remains the same for all models
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
    
    outlist[[i]]$stopoverfit <- temp
    outlist[[i]]$year <- yr
  }
  return(outlist)
}




####
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

# if modeled by gdd
counts <- counts[which(counts$cumdegday %% 50 == 0), ]
counts$Week <- counts$cumdegday / 50 + 1


# if modeled by Ordinal
counts$SiteDate <- strptime(paste(counts$Year, counts$Ordinal, sep = " "), format = "%Y %j")
counts$SiteDate <- as.Date(counts$SiteDate)
counts$Week <- week(counts$SiteDate)
counts$Week <- counts$Week - min(counts$Week) + 1
#overachieving volunteers going out more than once a week!
#choose first one, not averaging, easier to match with covariates
counts_uniq <- counts %>% 
  group_by(SiteID, Week) %>%
  arrange(Week) %>%
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

# dummy p covariate to make model work, but really just want common Np index
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


########################################################
# p.m <- "cov"
p.m <- "common"
w.m <- "cov"
mu.m <- "cov"
sigma.m <- "het"
# phi.m <- "logit.a"
phi.m <- "const"

Tries <- 3
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

