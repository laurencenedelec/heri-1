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

n_Estimation <- 2
n_MAX <- 500

##########
## Compare variation of importance of delta dom and epi in a specific setup
##########
param.list[[1]] <- list(n = rep(500, times = n_Estimation),
                       s = rep(500, times = n_Estimation),
                       u = rep(20, times = n_Estimation),
                       b = c(2,2),
                       delta_add = rep(1, times = n_Estimation),
                       delta_dom = seq(from = 0, to = 100, by = 100 / n_Estimation),
                       delta_epi = seq(from = 0, to = 100, by = 100 / n_Estimation), 
                       snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
                       variable = c("delta_dom", "delta_epi"))                   

param.list[[2]] <- list(n = rep(500, times = n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = seq(from = 0, to = 100, by = 100 / n_Estimation),
                        delta_epi = seq(from = 0, to = 100, by = 100 / n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
                        variable = c("delta_dom", "delta_epi"))    

param.list[[3]] <- list(n = rep(500, times = n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(1, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = seq(from = 0, to = 100, by = 100 / n_Estimation),
                        delta_epi = seq(from = 0, to = 100, by = 100 / n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
                        variable = c("delta_dom", "delta_epi"))    

##########
## When varying number of SnPS, also vary importance of dom and epi effect
##########
param.list[[4]] <- list(n = rep(500, times = n_Estimation),
                       s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                       u = rep(1, times = n_Estimation),
                       b = c(2,2),
                       delta_add = rep(1, times = n_Estimation),
                       delta_dom = rep(0, times = n_Estimation),
                       delta_epi = rep(0, times = n_Estimation),   
                       snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                        
                       variable = "s")

param.list[[5]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(1, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(1, times = n_Estimation),
                        delta_epi = rep(1, times = n_Estimation),    
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                       
                        variable = "s")

param.list[[6]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(1, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(10, times = n_Estimation),
                        delta_epi = rep(10, times = n_Estimation),    
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                       
                        variable = "s")

param.list[[7]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(1, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(100, times = n_Estimation),
                        delta_epi = rep(100, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                        
                        variable = "s")

##########
## When varying number of SnPS, also vary importance of dom and epi effect
##########
param.list[[8]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(0, times = n_Estimation),
                        delta_epi = rep(0, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                        
                        variable = "s")

param.list[[9]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(1, times = n_Estimation),
                        delta_epi = rep(1, times = n_Estimation),    
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                       
                        variable = "s")

param.list[[10]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(10, times = n_Estimation),
                        delta_epi = rep(10, times = n_Estimation),     
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                      
                        variable = "s")

param.list[[11]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(100, times = n_Estimation),
                        delta_epi = rep(100, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),                        
                        variable = "s")

##########
## When varying number of SnPS with low n, also vary importance of dom and epi effect
##########
param.list[[12]] <- list(n = rep(20, times = n_Estimation),
                        s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(0, times = n_Estimation),
                        delta_epi = rep(0, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),     
                        variable = "s")
 
param.list[[13]] <- list(n = rep(20, times = n_Estimation),
                         s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                         u = rep(10, times = n_Estimation),
                         b = c(2,2),
                         delta_add = rep(1, times = n_Estimation),
                         delta_dom = rep(1, times = n_Estimation),
                         delta_epi = rep(1, times = n_Estimation),   
                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),     
                         variable = "s")

param.list[[14]] <- list(n = rep(20, times = n_Estimation),
                         s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                         u = rep(10, times = n_Estimation),
                         b = c(2,2),
                         delta_add = rep(1, times = n_Estimation),
                         delta_dom = rep(10, times = n_Estimation),
                         delta_epi = rep(10, times = n_Estimation),   
                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),     
                         variable = "s")

param.list[[15]] <- list(n = rep(20, times = n_Estimation),
                         s = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                         u = rep(10, times = n_Estimation),
                         b = c(2,2),
                         delta_add = rep(1, times = n_Estimation),
                         delta_dom = rep(100, times = n_Estimation),
                         delta_epi = rep(100, times = n_Estimation),    
                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),    
                         variable = "s")

##########
## When varying number of real_coeff, also vary importance of dom and epi effect
##########
param.list[[16]] <- list(n = rep(n_MAX, times = n_Estimation),
                        s = rep(n_MAX, times = n_Estimation),
                        u = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(0, times = n_Estimation),
                        delta_epi = rep(0, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "u")
 
param.list[[17]] <- list(n = rep(n_MAX, times = n_Estimation),
                        s = rep(n_MAX, times = n_Estimation),
                        u = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(1, times = n_Estimation),
                        delta_epi = rep(1, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "u")

param.list[[18]] <- list(n = rep(n_MAX, times = n_Estimation),
                        s = rep(n_MAX, times = n_Estimation),
                        u = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(10, times = n_Estimation),
                        delta_epi = rep(10, times = n_Estimation),    
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T), 
                        variable = "u")

param.list[[19]] <- list(n = rep(n_MAX, times = n_Estimation),
                        s = rep(n_MAX, times = n_Estimation),
                        u = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(100, times = n_Estimation),
                        delta_epi = rep(100, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "u")

##########
## When varying number n, also vary importance of dom and epi effect
##########
param.list[[20]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(0, times = n_Estimation),
                        delta_epi = rep(0, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "n")

param.list[[21]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(1, times = n_Estimation),
                        delta_epi = rep(1, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "n")

param.list[[22]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(10, times = n_Estimation),
                        delta_epi = rep(10, times = n_Estimation),    
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T), 
                        variable = "n")

param.list[[23]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
                        s = rep(20, times = n_Estimation),
                        u = rep(10, times = n_Estimation),
                        b = c(2,2),
                        delta_add = rep(1, times = n_Estimation),
                        delta_dom = rep(100, times = n_Estimation),
                        delta_epi = rep(100, times = n_Estimation),   
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),  
                        variable = "n")


res <- list()
# Do estimation and generate parameters
for(i in 1:length(param.list)) {
    
    cat(paste0("= Param list: ", i, "/", length(param.list), "\n"))
    
    res[[i]] <- compare_dcor(n = param.list[[i]]$n, 
                             s = param.list[[i]]$s,
                             u = param.list[[i]]$u, 
                             b = param.list[[i]]$b, 
                             delta_add = param.list[[i]]$delta_add,
                             delta_dom = param.list[[i]]$delta_dom,
                             delta_epi = param.list[[i]]$delta_epi,
                             snps_value = param.list[[i]]$snps_value,
                             variable = param.list[[i]]$variable)
    
    cat(paste0("= End of ", i, "/", length(param.list), "\n"))
                       
}
