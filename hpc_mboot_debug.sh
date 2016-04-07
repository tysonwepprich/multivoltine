#!/bin/tcsh
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -W 15
setenv OMP_NUM_THREADS 4
source /usr/local/apps/R/R-2.15-parallel.csh
R CMD BATCH --vanilla  hpc_mbootstrap.R 
#BSUB -o out.%J
#BSUB -e err.%J