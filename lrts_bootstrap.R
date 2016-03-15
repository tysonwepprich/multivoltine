# Read data and do slurmGeneration
source('bootstrapMfunctions.R')

simdatafiles <- list.files("simDataGenMode/")
i <- 19 # choose species

SampleList <- readRDS(paste("simDataGenMode/", simdatafiles[i], sep = ""))
# data_file Rdata
dataIN <- c("SampleList")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = sample(seq(1:length(SampleList))))

# calculate null hypotheses for same species, different years
ZabSkip <- slurm_apply(f = SlurmGeneration, params = paramIN, 
                        cpus_per_node = 8, nodes = 4, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")


