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

#setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

### Settings
options(digits.secs=10)
Project <- SetupProject()


####################
### Comparison DCOR
####################
param.list <- list()

N_Estimation <- 50
N_MAX <- 500
param.list[[1]] <- list(N = rep(500, times = N_Estimation),
                       N_SNPS = rep(500, times = N_Estimation),
                       N_real_coeff = rep(20, times = N_Estimation),
                       b = c(2,2),
                       delta_add = rep(1, times = N_Estimation),
                       delta_dom = seq(from = 0, to = 1, by = 1 / N_Estimation),
                       delta_epi = seq(from = 0, to = 1, by = 1 / N_Estimation),   
                       variable = c("delta_dom", "delta_epi"))                   

param.list[[1]] <- list(N = rep(500, times = N_Estimation),
                       N_SNPS = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                       N_real_coeff = rep(1, times = N_Estimation),
                       b = c(2,2),
                       delta_add = rep(1, times = N_Estimation),
                       delta_dom = rep(1, times = N_Estimation),
                       delta_epi = rep(1, times = N_Estimation),                        
                       variable = "N_SNPS")
 
param.list[[2]] <- list(N = rep(500, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        N_real_coeff = rep(10, times = N_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = N_Estimation),
                        delta_dom = rep(1, times = N_Estimation),
                        delta_epi = rep(1, times = N_Estimation),   
                        variable = "N_SNPS")
 
param.list[[3]] <- list(N = rep(20, times = N_Estimation),
                        N_SNPS = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        N_real_coeff = rep(10, times = N_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = N_Estimation),
                        delta_dom = rep(1, times = N_Estimation),
                        delta_epi = rep(1, times = N_Estimation),     
                        variable = "N_SNPS")
 
param.list[[4]] <- list(N = rep(N_MAX, times = N_Estimation),
                        N_SNPS = rep(N_MAX, times = N_Estimation),
                        N_real_coeff = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = N_Estimation),
                        delta_dom = rep(1, times = N_Estimation),
                        delta_epi = rep(1, times = N_Estimation),  
                        variable = "N_real_coeff")
 
param.list[[5]] <- list(N = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        N_SNPS = rep(20, times = N_Estimation),
                        N_real_coeff = rep(10, times = N_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = N_Estimation),
                        delta_dom = rep(1, times = N_Estimation),
                        delta_epi = rep(1, times = N_Estimation),  
                        variable = "N")

param.list[[6]] <- list(N = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        N_SNPS = seq(from = 10, to = N_MAX, by = N_MAX/N_Estimation),
                        N_real_coeff = rep(10, times = N_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = N_Estimation),
                        delta_dom = rep(1, times = N_Estimation),
                        delta_epi = rep(1, times = N_Estimation),  
                        variable = c("N", "N_SNPS"))

res <- list()
# Do estimation and generate parameters
for(i in 1:length(param.list)) {
    
    cat(paste0("= Param list: ", i, "/", length(param.list), "\n"))
    
    res[[i]] <- compare_dcor(N = param.list[[i]]$N, 
                             N_SNPS = param.list[[i]]$N_SNPS,
                             N_real_coeff = param.list[[i]]$N_real_coeff, 
                             b = param.list[[i]]$b, 
                             delta_add = param.list[[i]]$delta_add,
                             delta_dom = param.list[[i]]$delta_dom,
                             delta_epi = param.list[[i]]$delta_epi,
                             variable = param.list[[i]]$variable)
    
    cat(paste0("= End of ", i, "/", length(param.list), "\n"))
                       
}
