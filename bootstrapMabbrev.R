# abbreviated version of bootstrampM.R
# run one M null hypothesis at a time, to prevent unnecessary computation

source('bootstrapMfunctions.R')


# load data
count_array <- readRDS('count_array_expanded.rds')
cov_array <- readRDS('covariates_array_expanded.rds')
cov_sites <- readRDS('covariates_sites_expanded.rds')
species_list <- readRDS('species_expanded.rds')
# 
# count_array <- readRDS('count_array.2009.rds')
# cov_array <- readRDS('covariates_array.2009.rds')
# cov_sites <- readRDS('covariates_sites.2009.rds')
# species_list_2009 <- readRDS('species.2009.rds')

# choose parameter ranges
species <- c(6, 20, 21, 23, 25, 27, 30, 31, 32, 33, 34) # corresponds to row in species_list
raw_cutoff <- 10 
p_cov1 <- 7 # Select detection covariates here (1:7 possible)
p_cov2 <- "none" # c("none", 1:6) # Select detection covariates here (1:7 possible)
site_covs <- "AnnGDD" # c("common", "AnnGDD", "SprGDD", "lat") # for mu, w 
M <- c(1, 1, 1, 2, 3) #number of broods to estimate
sigma.m <- "het" #  c("het", "hom")
phi.m <- "const" # c("const", "logit.a")

params <- expand.grid(species, raw_cutoff, p_cov1, p_cov2, 
                      site_covs, M, sigma.m, phi.m,
                      stringsAsFactors = FALSE)
names(params) <- c("species", "raw_cutoff", "p_cov1", "p_cov2", 
                   "site_covs", "M", "sigma.m", "phi.m")

# data_file Rdata
dataIN <- c("count_array", "cov_array", "cov_sites", "params")
save(list = dataIN, file = "dataIN.RData")

# simple param file for slurm.apply
paramIN <- data.frame(nRun = seq(1:nrow(params)))

# calculate null hypotheses for M for different species
baseline <- slurm_apply(f = SlurmCovs, params = paramIN, 
                        cpus_per_node = 8, nodes = 4, 
                        data_file = "dataIN.RData", 
                        # pkgs = c("devtools", "msm", "rslurm", "StopoverCode"), 
                        output = "raw")





