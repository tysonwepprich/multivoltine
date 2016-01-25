# Hopefully final draft of stopover model analysis

# 1. By species, constrain M to best values (by manually looking at phenograms)
# 2. Auto-bootstrap M selection for base M and M+1 for each year
# 3. Run bootstrap to get parameter estimates for best M for each year

# Difficult thing is optimizing for slurm. May not be able to have one big function,
# since some steps may need many cores, some few

# file needs to be edited depending on whether Windows or Linux for package loading
source('bootstrapMfunctions.R')

allSpecies <- read.csv("data/MultivoltineSpecies.csv", header = TRUE)

# choose your species
i <- 19
species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i]

# somewhat unwieldly, list with each year as a list of 4 (year, counts, surv_covs, site_covs)
dat <- SpeciesData(species)

# for each species, select parameters
# how much data available for modeling?

count_cutoff <- 5
surv_cutoff <- 3
data_avail <- data.frame()
for (j in 1:length(dat)){
  list_index <- j
  temp <- dat[[j]]
  year <- temp[[1]]
  counts <- temp[[2]]
  sum_count <- apply(counts, 1, sum, na.rm = TRUE)
  surv_present <- apply(counts, 1, function(x) length(which(x > 0)))
  cutoff_met <- length(which(sum_count >= count_cutoff))
  survs_met <- length(which(surv_present >= surv_cutoff))
  both <- cbind(sum_count, surv_present)
  both_met <- length(which(apply(both, 1, function(x, count_cutoff, surv_cutoff) x[1] >= count_cutoff & 
                                   x[2] >= surv_cutoff, count_cutoff = count_cutoff, surv_cutoff = surv_cutoff)))
  new_row <- data.frame(count_cutoff, surv_cutoff, year, list_index, cutoff_met, survs_met, both_met)
  data_avail <- rbind(data_avail, new_row)
}

list_index_min_data <- unique(data_avail$list_index[data_avail$both_met >= 10])

# choose parameter ranges
raw_cutoff <- 5 # c(5, 10)
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- c("AnnGDD", "lat") # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(minBrood:maxBrood) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, list_index_min_data, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "list_index", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("dat", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

cl <- makeCluster(7)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmCovs)
stopCluster(cl)

saveRDS(test, file = "ZabSkipCovs.rds")

# choose your species
i <- 18
species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i]

# somewhat unwieldly, list with each year as a list of 4 (year, counts, surv_covs, site_covs)
dat <- SpeciesData(species)

# for each species, select parameters
# how much data available for modeling?

count_cutoff <- 5
surv_cutoff <- 3
data_avail <- data.frame()
for (j in 1:length(dat)){
  list_index <- j
  temp <- dat[[j]]
  year <- temp[[1]]
  counts <- temp[[2]]
  sum_count <- apply(counts, 1, sum, na.rm = TRUE)
  surv_present <- apply(counts, 1, function(x) length(which(x > 0)))
  cutoff_met <- length(which(sum_count >= count_cutoff))
  survs_met <- length(which(surv_present >= surv_cutoff))
  both <- cbind(sum_count, surv_present)
  both_met <- length(which(apply(both, 1, function(x, count_cutoff, surv_cutoff) x[1] >= count_cutoff & 
                                   x[2] >= surv_cutoff, count_cutoff = count_cutoff, surv_cutoff = surv_cutoff)))
  new_row <- data.frame(count_cutoff, surv_cutoff, year, list_index, cutoff_met, survs_met, both_met)
  data_avail <- rbind(data_avail, new_row)
}

list_index_min_data <- unique(data_avail$list_index[data_avail$both_met >= 10])

# choose parameter ranges
raw_cutoff <- 5 # c(5, 10)
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- c("AnnGDD", "lat") # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(minBrood:maxBrood) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, list_index_min_data, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "list_index", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("dat", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

cl <- makeCluster(7)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmCovs)
stopCluster(cl)

saveRDS(test, file = "WIDWcov.rds")

# choose your species
i <- 17
species <- allSpecies$CommonName[i]
minBrood <- allSpecies$MinBrood[i]
maxBrood <- allSpecies$MaxBrood[i]

# somewhat unwieldly, list with each year as a list of 4 (year, counts, surv_covs, site_covs)
dat <- SpeciesData(species)

# for each species, select parameters
# how much data available for modeling?

count_cutoff <- 5
surv_cutoff <- 3
data_avail <- data.frame()
for (j in 1:length(dat)){
  list_index <- j
  temp <- dat[[j]]
  year <- temp[[1]]
  counts <- temp[[2]]
  sum_count <- apply(counts, 1, sum, na.rm = TRUE)
  surv_present <- apply(counts, 1, function(x) length(which(x > 0)))
  cutoff_met <- length(which(sum_count >= count_cutoff))
  survs_met <- length(which(surv_present >= surv_cutoff))
  both <- cbind(sum_count, surv_present)
  both_met <- length(which(apply(both, 1, function(x, count_cutoff, surv_cutoff) x[1] >= count_cutoff & 
                                   x[2] >= surv_cutoff, count_cutoff = count_cutoff, surv_cutoff = surv_cutoff)))
  new_row <- data.frame(count_cutoff, surv_cutoff, year, list_index, cutoff_met, survs_met, both_met)
  data_avail <- rbind(data_avail, new_row)
}

list_index_min_data <- unique(data_avail$list_index[data_avail$both_met >= 10])

# choose parameter ranges
raw_cutoff <- 5 # c(5, 10)
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- c("AnnGDD", "lat") # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(minBrood:maxBrood) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, list_index_min_data, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "list_index", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("dat", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

cl <- makeCluster(7)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmCovs)
stopCluster(cl)

saveRDS(test, file = "ViceroyCov.rds")
