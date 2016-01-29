#parallel benchmarking
source('bootstrapMfunctions.R')



SampleList <- readRDS("simDataGenMode/NorthernBroken-Dashsimdata.rds")
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

strat <- c("parLapply", "parLapplyLB","parLapplyLB2", "rev_parLapply", "rev_parLapplyLB","rev_parLapplyLB2")
speedtest <- data.frame(strat, time = NA)

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(from = 1, to =length(SampleList)))


# parLapply
speedtest$time[1] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})

saveRDS(test, file = "NBD_bootstrapM.rds")

# parLapplyLB
speedtest$time[2] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapplyLB(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})

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
speedtest$time[3] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapplyLB2(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})

# reverse order, complex models and longer times at front
paramIN <- data.frame(nRun = seq(from = length(SampleList), to = 1))


# rev_parLapply
speedtest$time[4] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapply(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})

# rev_parLapplyLB
speedtest$time[5] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapplyLB(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})


# rev_parLapplyLB2
speedtest$time[6] <- system.time({cl <- makeCluster(16)
clusterEvalQ(cl, {
  library(devtools)
  library(msm)
  library(dplyr)
  # library(StopoverCode) #on linux
  devtools::load_all("StopoverCode", recompile = TRUE) # on windows
  load("dataIN.RData")
})
test <- parLapplyLB2(cl, paramIN$nRun, SlurmGeneration)
stopCluster(cl)})

saveRDS(speedtest, file = "speedtest.rds")