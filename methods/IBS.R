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
####  IBS
####
####**********************************************************************

# library used in this R script
#library()

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
get_estimate_IBS <- function(V, phi) {

    h <- estimate_heritability(V, as.vector(phi))
    
    list(h = h$heritability)
  
}

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
get_estimate_IBS_dcov <- function(V, phi) {
    
    h <- estimate_heritability_dcov(V, phi)
    
    list(h = h$heritability)
    
}


#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
get_estimate_IBS_dcov_LN <- function(V, phi) {
    
    h <- estimate_heritability_dcov_LN(V, phi)
    
    list(h = h$heritability)
    
}