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
####  Heritability estimate through linear regression
####
####**********************************************************************

# library used in this R script
#library()

#' Estimate heritability
#'
#' @title Estimate heritability
#' @param V matrix built from P,K
#' @param phi 2*Kinship
#' @return heritability + p.value obtained from lm
#' @author Julien Duvanel
#' @export
estimate_heritability <- function(V, phi) {
  
  # Get data for the specific phenotype
  Y <- V
  Y.length <- length(Y)
  
  Y.square <- Y %*% t(Y)
  
  # We need this as a vector
  Y.square.vec <- as.vector(Y.square)
  
  # Compute the heritability using a linear regression
  lm.heritability <- lm(Y.square.vec ~ phi)
  
  heritability <- summary(lm.heritability)$coef[2,1]
  heritability.p.value <- summary(lm.heritability)$coef[2,4]
  
  # Return
  list(heritability = heritability, p.value = heritability.p.value)
  
}