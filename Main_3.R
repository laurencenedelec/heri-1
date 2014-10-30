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

# Parallel computing
suppressPackageStartupMessages(library(doMC))
registerDoMC(4)

### Clean workspace, set directory and load functions
rm(list = ls())

setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
#setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

### Settings
options(digits.secs=10)
Project <- SetupProject()

N <- c(100, 1000)
N_SNPS <- c(100, 100)
N_real_coeff <- c(20, 20)
b <- c(2,2)

compare_dcor(N, N_SNPS,N_real_coeff, build_SNPs_matrix, b, c("N", "N_SNPs"))
