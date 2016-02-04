#parallel benchmarking
source('bootstrapMfunctions.R')

ncores <- detectCores()

SampleList <- readRDS("simDataGenMode/NorthernBroken-Dashsimdata.rds")
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

strat <- c("parLapply", "parLapplyLB","parLapplyLB2", "rev_parLapply", "rev_parLapplyLB","rev_parLapplyLB2")
speedtest <- data.frame(strat, time = NA)

# split param file for 6 tests

a <- seq(from = 1, to = length(SampleList), by = 1)

paramIN1 <- data.frame(nRun = a[seq(1, length(a), 6)])
paramIN2 <- data.frame(nRun = a[seq(2, length(a), 6)])
paramIN3 <- data.frame(nRun = a[seq(3, length(a), 6)])
paramIN4 <- data.frame(nRun = sort( a[seq(4, length(a), 6)], decreasing = TRUE))
paramIN5 <- data.frame(nRun = sort( a[seq(5, length(a), 6)], decreasing = TRUE))
paramIN6 <- data.frame(nRun = sort( a[seq(6, length(a), 6)], decreasing = TRUE))


# parLapply
speedtest$time[1] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  library(StopoverCode) #on linux
 # devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test1 <- parLapply(cl, paramIN1$nRun, SlurmGeneration)
stopCluster(cl)})

saveRDS(test, file = "NBD_bootstrapM1.rds")

# parLapplyLB
speedtest$time[2] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
   library(StopoverCode) #on linux
  #devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test2 <- parLapplyLB(cl, paramIN2$nRun, SlurmGeneration)
stopCluster(cl)})
saveRDS(test2, file = "NBD_bootstrapM2.rds")

# alternative parLapply load balancing function
parLapplyLB2 <- function(cl, x, fun, ...)
{
  LB.init <- function(fun, ...)
  {
    assign(".LB.fun", fun, pos=globalenv())
    assign(".LB.args", list(...), pos=globalenv())
    NULL
  }
  
  LB.worker <- function(x) do.call(".LB.fun", c(list(x), .LB.args))
  
  clusterCall(cl, LB.init, fun, ...)
  r <- clusterApplyLB(cl, x, LB.worker)
  clusterEvalQ(cl, rm(".LB.fun", ".LB.args", pos=globalenv()))
  r
} 

# parLapplyLB2
speedtest$time[3] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
   library(StopoverCode) #on linux
  #devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test3 <- parLapplyLB2(cl, paramIN3$nRun, SlurmGeneration)
stopCluster(cl)})
saveRDS(test3, file = "NBD_bootstrapM3.rds")


# rev_parLapply
speedtest$time[4] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
   library(StopoverCode) #on linux
  #devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test4 <- parLapply(cl, paramIN4$nRun, SlurmGeneration)
stopCluster(cl)})
saveRDS(test4, file = "NBD_bootstrapM4.rds")


# rev_parLapplyLB
speedtest$time[5] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
   library(StopoverCode) #on linux
  #devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test5 <- parLapplyLB(cl, paramIN5$nRun, SlurmGeneration)
stopCluster(cl)})
saveRDS(test5, file = "NBD_bootstrapM5.rds")


# rev_parLapplyLB2
speedtest$time[6] <- system.time({cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
   library(StopoverCode) #on linux
  #devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test6 <- parLapplyLB2(cl, paramIN6$nRun, SlurmGeneration)
stopCluster(cl)})
saveRDS(test6, file = "NBD_bootstrapM6.rds")


saveRDS(speedtest, file = "speedtest.rds")

