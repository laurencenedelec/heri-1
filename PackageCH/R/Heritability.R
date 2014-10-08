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

#' Estimate heritability
#'
#' @title Estimate heritability
#' @param V matrix built from P,K
#' @param phi 2*Kinship
#' @return heritability + p.value obtained from lm
#' @author Julien Duvanel
#' @export
estimate_heritability_dcov <- function(V, phi) {
    
    # Get data for the specific phenotype
    
    Y.distance.carre <- matrix(rep(V^2,nrow(V)),
                               ncol = nrow(V)) + 
        t(matrix(rep(V^2,nrow(V)),
                 ncol = nrow(V))) - 
        2 * V %*%t(V)
    
    Y.distance <- Y.distance.carre^(1/2)
    
    phi.matrix <- (2 * phi)^(1/2)

    heritability <- dcov(Y.distance, phi.matrix) / dcov(phi.matrix, phi.matrix)
    
    # Return
    list(heritability = heritability)
    
}