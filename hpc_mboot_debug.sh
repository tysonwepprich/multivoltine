#!/bin/tcsh
#BSUB -n 4
#BUSB -W 15
source /usr/local/apps/R/R-2.15-parallel.csh
R CMD BATCH --vanilla  hpc_mbootstrap.r 
#BSUB -o out.%J
#BSUB -e err.%J