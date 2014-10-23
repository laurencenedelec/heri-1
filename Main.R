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
registerDoMC(1)

### Clean workspace, set directory and load functions
rm(list = ls())

setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
#setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

# Load every methods and data
source("methods/GCTA.R")
source("methods/loadGCTA.R")

source("methods/TheoKin.R")
source("methods/loadTheoKin.R")

source("methods/IBS.R")
source("methods/loadIBS.R")

source("methods/IBD.R")
source("methods/loadIBD.R")

### Settings
options(digits.secs=10)
Project <- SetupProject()

### List of parameters

# Load phenotypes
# The format is strict:
# two first columns must be: famid, id
# then phenotypes (and we give phenotypes id corresponding to the column below)
# source("methods/loadPhenotypes.R")
# 
# param.list <- list()
# param.list[[1]] <- list(method.to.test = c(
#                                            "GCTA",
#                                            "TheoKin",
#                                            "IBS",
#                                            "IBD"
#                                            ),
#                         model.to.test = c(
#                                             "",
#                                             "_dcov",
#                                             "_dcov_LN"
#                                         ),
#                         phenotypes.id = c(3:5)
#                        )
# 
# try(DoBatchEstimation(param.list = param.list))

# Deuxieme estimation (avec les vrais traits)
source("methods/loadPhenotypes_real.R")

param.list <- list()
param.list[[1]] <- list(method.to.test = c(
                                           "GCTA"#,
                                           #"TheoKin",
                                           #"IBS",
                                           #"IBD"
                                           ),
                        model.to.test = c(
                                            #"", 
                                            #"_dcov",
                                            #"_dcov_LN"
                                            "_PlotSimilarity"
                                        ),
                        phenotypes.id = c(4:5)
                       )

try(DoBatchEstimation(param.list = param.list))
