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
####  Theoretical Kinship
####
####**********************************************************************

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
get_estimate_TheoKin <- function(V, phi) {

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
get_estimate_TheoKin_dcov <- function(V, phi) {
    
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
get_estimate_TheoKin_dcov_LN <- function(V, phi) {
    
    h <- estimate_heritability_dcov_LN(V, phi)
    
    list(h = h$heritability)
    
}