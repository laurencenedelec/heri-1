####**********************************************************************
####  Written and Developed by: 
####**********************************************************************
####
####    Julien Duvanel, Copyright 2014
####    email: duvanel@stanford.edu
####
####**********************************************************************
####**********************************************************************

####**********************************************************************
####
####  Main file
####
####**********************************************************************

suppressPackageStartupMessages(library(devtools))

### Clean workspace, set directory and load functions
rm(list = ls())

setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
#setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

### Settings
options(digits.secs=10)
Project <- SetupProject()

N_Estimation <- 100
N <- rep(400, times = N_Estimation)
N_SNPS <- seq(from = 10, to = 1000, by = (1000-10+10)/N_Estimation)
N_real_coeff <- rep(0, times = N_Estimation)
b <- c(2,2)
variable <- "N_SNPS"

res <- compare_dcor(N, N_SNPS,N_real_coeff, build_SNPs_matrix, b, variable)
