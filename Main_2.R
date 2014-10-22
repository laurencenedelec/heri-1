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

#setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

### Settings
options(digits.secs=10)
Project <- SetupProject()

# Used to change name when we save file
datetime.stamp <- format(Sys.time(), "%d%m%Y_%H%M%S")

# Load matrix G
G <- build_matrix_G_manually()

# Load matrix P
load("data/Phenotype_Simulated_14102014_145630.RData")

# Estimate heritability
h <- HeritabilityEstimation(G, P.raw$fake.trait)

save(list = c("G", "h"), file = paste0("results.", datetime.stamp, ".RData"))
