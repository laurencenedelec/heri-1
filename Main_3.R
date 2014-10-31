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


####################
### Comparison DCOR
####################

param.list <- list()

N_Estimation <- 100
param.list[[1]] <- list(N = rep(400, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = 1000, by = (1000-10+10)/N_Estimation),
                        N_real_coeff = rep(0, times = N_Estimation),
                        b = c(2,2),
                        variable = "N_SNPS",
                        AddNoise = TRUE)
                        
param.list[[2]] <- list(N = rep(400, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = 1000, by = (1000-10+10)/N_Estimation),
                        N_real_coeff = rep(0, times = N_Estimation),
                        b = c(2,2),
                        variable = "N_SNPS",
                        AddNoise = FALSE)                        

param.list[[3]] <- list(N = rep(400, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = 1000, by = (1000-10+10)/N_Estimation),
                        N_real_coeff = rep(1, times = N_Estimation),
                        b = c(2,2),
                        variable = "N_SNPS",
                        AddNoise = TRUE)

param.list[[4]] <- list(N = rep(400, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = 1000, by = (1000-10+10)/N_Estimation),
                        N_real_coeff = rep(1, times = N_Estimation),
                        b = c(2,2),
                        variable = "N_SNPS",
                        AddNoise = FALSE)

param.list[[1]] <- list(N = rep(400, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = 1000, by = 1000/N_Estimation),
                        N_real_coeff = rep(10, times = N_Estimation),
                        b = c(2,2),
                        variable = "N_SNPS",
                        AddNoise = TRUE)



res <- list()
# Do estimation and generate parameters
for(i in 1:length(param.list)) {
    
    cat(paste0("= Param list: ", i, "/", length(param.list), "\n"))
    
    res[[i]] <- compare_dcor(N = param.list[[i]]$N, 
                        N_SNPS = param.list[[i]]$N_SNPS,
                        N_real_coeff = param.list[[i]]$N_real_coeff, 
                        get_snps_matrix = build_SNPs_matrix, 
                        get_alpha = build_alpha,
                        b = param.list[[i]]$b, 
                        variable = param.list[[i]]$variable,
                        AddNoise = param.list[[i]]$AddNoise)
        
    cat(paste0("= End of ", i, "/", length(param.list), "\n"))
                       
}
