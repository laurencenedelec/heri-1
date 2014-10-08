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
####  GCTA
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
get_estimate_GCTA <- function(V, phi) {

    h <- estimate_heritability(V, phi)

    list(h = h$heritability)
  
}

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
get_estimate_GCTA_dcov <- function(V, phi) {
    
    h <- estimate_heritability_dcov(V, phi)

    list(h = h$heritability)
    
}