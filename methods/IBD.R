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
####  IBD
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
get_estimate_IBD <- function(P, V, phi, dcov = F) {

  if(dcov == TRUE) {
    h <- estimate_heritability_dcov(V, phi)
  } else {
    h <- estimate_heritability(V, phi)
  }
  list(h = h$heritability)
  
}
