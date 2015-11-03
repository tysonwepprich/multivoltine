# finding results when you clear the workspace and get_slurm_out no longer works
slurm_codes <- c("slr3198", "slr9298", "slr448", "slr4166", "slr9994")
slurm_out <- list()

for (j in 1:length(slurm_codes)){
  missing_files <- c()
  tmpEnv <- new.env()
  for (i in 0:9) {
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

saveRDS(slurm_out, "modelResultsList.rds")

outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]]$pars
  out$model <- i
  out$ll.val <- outList[[i]]$ll.val
  if (is.na(out$ll.val)){
    out$npar <- NA
  }else{
    out$npar <- outList[[i]]$npar
  }
  out$time <- as.double(outList[[i]]$time, units = "mins")
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)


# forgot to label output with year, which is necessary since species labels change
# or need to use Common Name instead of index within count_array

outDF$year <- 2008
outDF$year[145:624] <- 2009

outDF2008 <- outDF[outDF$year == 2008, ]
outDF2008$species <- mapvalues(outDF2008$species, 
            from = c(9:21, 23:25, 27, 28, 30:32),
            to = species_list$CommonName[c(9:21, 23:25, 27, 28, 30:32)])
outDF2009 <- outDF[outDF$year == 2009, ]
outDF2009$species <- mapvalues(outDF2009$species, 
                               from = c(9:16, 18, 20),
                               to = species_list_2009$CommonName[c(9:16, 18, 20)])
outDF <- rbind(outDF2008, outDF2009)
saveRDS(outDF, "bestCovPrelimResults.rds")

spec <- unique(outDF$species)

covTest <- outDF[is.na(outDF$ll.val) == FALSE, ]
covTest$AIC <- -2 * covTest$ll.val + 2 * covTest$npar

test <- covTest %>%
  group_by(species, year) %>%
  mutate(weight1 = exp(-0.5 * (AIC - min(AIC)))) %>%
  mutate(weight = weight1 / sum(weight1))

test2008 <- test %>%
  select(-raw_cutoff, -weight1, -p_cov1, -p_cov2) %>%
  filter(species == spec[21] & year == 2008) %>%
  arrange(AIC) %>%
  data.frame()
test2009 <- test %>%
  select(-raw_cutoff, -weight1, -p_cov1, -p_cov2) %>%
  filter(species == spec[21] & year == 2009) %>%
  arrange(AIC) %>%
  data.frame()
head(test2008)
head(test2009)
