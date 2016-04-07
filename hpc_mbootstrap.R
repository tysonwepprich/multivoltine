# Bootstrap test for M (number of broods)
# Script to run on NCSU HPC
source('bootstrapMfunctions.R')


rdsfiles <- list.files("simDataGenMode/")
j <- 1

sp <- unlist(strsplit(rdsfiles[j], split = "cov", fixed = TRUE))[1]

SampleList <- readRDS(paste("simDataGenMode/", sp, sep = ""))
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = sample(seq(1:length(SampleList)))) # random nRun so split even for parallel


cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  library(StopoverCode) #on linux
  # devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
time <- system.time({test <- parLapply(cl, paramIN$nRun, SlurmGenerationP1)})
stopCluster(cl)

saveRDS(test, file = paste(sp, "BSmod.rds", sep = ""))
