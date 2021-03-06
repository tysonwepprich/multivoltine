# take slurm covariates output and put into data files by species 
# this data will then be ready for simulations for mixture model mode testing

# also should extract needed data from output from non-slurm SlurvCov analysis

allSpecies <- read.csv("data/MultivoltineSpecies.csv", header = TRUE)

# choose your species
# i <- 16 #slr4158 
# i <- 15 #slr7296
# i <- 14 #slr7389
# i <- 13 #slr1826
# i <- 12 #
# i <- 10 #
# i <- 2 
i <- 3 #slr2152, #slr2085 duplicated accidentally, slr1965, slr2023 also ETS??

########################
# extract data from SlurmCov results from parlapply (non-slurm)
# temp <- readRDS("slurmCovOutput/otherResults/ZabSkipCovs.rds")
temp <- readRDS("SSSKIP_covtest.rds")

test <- do.call(rbind, lapply(temp, function(x) length(x[[1]]))) # extra layer of list

outList <- temp
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]][[1]]$pars
  out$model <- i
  out$ll.val <- outList[[i]][[1]]$ll.val
  if (is.na(out$ll.val)){
    out$npar <- NA
  }else{
    out$npar <- outList[[i]][[1]]$npar
  }
  out$time <- as.double(outList[[i]][[1]]$time, units = "mins")
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
baselineDF <- outDF
#############################################


# extract data from SlurmCov results
slurm_codes <- c("slr3941")
slurm_out <- list()
# setwd("slurmCovOutput/sesyncResults")

for (j in 1:length(slurm_codes)){
  missing_files <- c()
  tmpEnv <- new.env()
  for (i in 0:11) {
    fname <- paste0(slurm_codes[j], "_", i, 
                    ".RData")
    if (fname %in% dir()) {
      load(fname, envir = tmpEnv)
      slurm_out <- c(slurm_out, get(".rslurm_result", 
                                    envir = tmpEnv))
    }
    else {
      missing_files <- c(missing_files, fname)
    }
  }
}
test <- do.call(rbind, lapply(slurm_out, function(x) length(x)))
# setwd("../../")

outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  if (length(outList[[i]]) == 1){
    out <- NA
  }else{
    out <- outList[[i]]$pars
    out$model <- i
    out$ll.val <- outList[[i]]$ll.val
    if (is.na(out$ll.val)){
      out$npar <- NA
    }else{
      out$npar <- outList[[i]]$npar
    }
    out$time <- as.double(outList[[i]]$time, units = "mins")
  }
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
baselineDF <- outDF
###############################################


# AnnGDD or latitude for site covariate? Judge with AIC

covTest <- baselineDF[is.na(baselineDF$ll.val) == FALSE, ]
covTest$AIC <- -2 * covTest$ll.val + 2 * covTest$npar
# covTest$AICc <- -2 * covTest$ll.val + 2 * covTest$npar
# minAIC <- min(covTest$AIC)
# covTest$weight1 <- exp(-0.5 * (covTest$AIC - minAIC))
# sumweight <- sum(covTest$weight1)
# covTest$weight <- covTest$weight1 / sumweight

test <- covTest %>%
  group_by(list_index, M, raw_cutoff) %>%
  mutate(weight1 = exp(-0.5 * (AIC - min(AIC)))) %>%
  mutate(weight = round(weight1 / sum(weight1), 3))

test2 <- test %>%
  select(-species, -weight1) %>%
  # filter(species == 11) %>%
  arrange(list_index, M, AIC) %>%
  data.frame()

# on average over years, does latitude or AnnGDD have higher weight 
test3 <- test2 %>%
  group_by(M, p_cov1, p_cov2, p_cov3, p_cov4, sigma.m) %>%
  summarise(mean_weight = mean(weight))

a <- test2[, c("list_index", "site_covs", "M", "weight")]
a$weight <- round(a$weight, 3)

# select results from baseline for simulations
# simulations/model fits to see how many modes in mixture model

# index <- which(baselineDF$M == 2 & baselineDF$site_covs == "lat") #Spicebush
# index <- which(baselineDF$M == 2 & baselineDF$site_covs == "lat") #Pecks
index <- which(baselineDF$M == 2 & baselineDF$site_covs == "lat") #ETS


slurm_out2 <- slurm_out[c(index)]

