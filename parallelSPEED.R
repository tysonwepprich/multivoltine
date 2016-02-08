#parallel benchmarking
source('bootstrapMfunctions.R')

ncores <- 8

SampleList <- readRDS("simDataGenMode/NorthernBroken-Dashsimdata.rds")
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

strat <- c("parLapply", "parLapplyLB","clusterApplyLB", "rev_parLapply", "rev_parLapplyLB","rev_clusterApplyLB")
speedtest <- as.list(strat)

# split param file for 6 tests

#for testing, just do 36 of SampleList

a <- seq(from = 1, to = length(SampleList), by = 1)

paramIN1 <- data.frame(nRun = a[seq(1, length(a), 6)])
paramIN2 <- data.frame(nRun = a[seq(2, length(a), 6)])
paramIN3 <- data.frame(nRun = a[seq(3, length(a), 6)])
paramIN4 <- data.frame(nRun = sort( a[seq(4, length(a), 6)], decreasing = TRUE))
paramIN5 <- data.frame(nRun = sort( a[seq(5, length(a), 6)], decreasing = TRUE))
paramIN6 <- data.frame(nRun = sort( a[seq(6, length(a), 6)], decreasing = TRUE))

cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})

# parLapply
speedtest[[1]][2] <- system.time({
test1 <- parLapply(cl, paramIN1$nRun, SlurmGeneration)
})
saveRDS(test1, file = "NBD_bootstrapM1.rds")

# parLapplyLB
speedtest$time[[2]][2] <- system.time({
test2 <- parLapplyLB(cl, paramIN2$nRun, SlurmGeneration)
})
saveRDS(test2, file = "NBD_bootstrapM2.rds")


# clusterApplyLB
speedtest$time[[3]][2] <- system.time({
test3 <- clusterApplyLB(cl, paramIN3$nRun, SlurmGeneration)
})
saveRDS(test3, file = "NBD_bootstrapM3.rds")


# rev_parLapply
speedtest$time[[4]][2] <- system.time({
test4 <- parLapply(cl, paramIN4$nRun, SlurmGeneration)
})
saveRDS(test4, file = "NBD_bootstrapM4.rds")


# rev_parLapplyLB
speedtest$time[[5]][2] <- system.time({
test5 <- parLapplyLB(cl, paramIN5$nRun, SlurmGeneration)
})
saveRDS(test5, file = "NBD_bootstrapM5.rds")


# rev_clusterApplyLB
speedtest$time[[6]][2] <- system.time({
test6 <- clusterApplyLB(cl, paramIN6$nRun, SlurmGeneration)
})
saveRDS(test6, file = "NBD_bootstrapM6.rds")

stopCluster(cl)

saveRDS(speedtest, file = "speedtest.rds")

