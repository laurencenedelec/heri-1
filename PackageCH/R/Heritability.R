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
library(energy)
library(ggplot2)

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
    
    Y.distance <- Y.distance.carre
    
    phi.matrix <- phi

    heritability <- dcov(Y.distance, phi.matrix) / dcov(phi.matrix, phi.matrix)
    
    # Return
    list(heritability = heritability)
    
}

#' Estimate heritability with dcov / LN version
#'
#' @title Estimate heritability with dcov / LN version
#' @param V matrix built from P,K
#' @param phi 2*Kinship
#' @return heritability + p.value obtained from lm
#' @author Julien Duvanel
#' @export
estimate_heritability_dcov_LN <- function(V, phi) {
    
    N <- length(V)
    
    # Get data for the specific phenotype   
    Y.distance.carre <- matrix(rep(V^2,nrow(V)),
                               ncol = nrow(V)) + 
        t(matrix(rep(V^2,nrow(V)),
                 ncol = nrow(V))) - 
        2 * V %*%t(V)
    
    Y.distance <- sqrt(Y.distance.carre)
    
    phi.matrix <- sqrt(N*(2-2*phi))
    #phi.matrix <- phi

    heritability <- dcov(Y.distance, phi.matrix) / (sqrt(dcov(phi.matrix, phi.matrix)) * sqrt(dcov(Y.distance, Y.distance)))
    
    # Return
    list(heritability = heritability)
    
}

#' Estimate heritability with dcov / LN version
#'
#' @title Estimate heritability with dcov / LN version
#' @param V matrix built from P,K
#' @param phi 2*Kinship
#' @return heritability + p.value obtained from lm
#' @author Julien Duvanel
#' @export
estimate_heritability_PlotSimilarity <- function(V, phi) {
    
    # datetime stamp (to save files)
    datetime.stamp <- format(Sys.time(), 
                             "%d%m%Y_%H%M%S")
    
    # Get data for the specific phenotype
    Y <- V
    Y.square <- Y %*% t(Y)
    
    # We need this as a vector
    Y.square.vec <- as.vector(Y.square)

    # heri
    lm.heritability <- lm(Y.square.vec ~ as.vector(phi))
    heritability <- summary(lm.heritability)$coef[2,1]
    
    df <- data.frame(phenotype.similarity = Y.square.vec,
                     genotype.similarity = as.vector(phi))
    
    p <- ggplot(data = df[sample(1:nrow(df), 0.05 * nrow(df), replace=FALSE),],
                aes(y = phenotype.similarity, 
                    x = genotype.similarity)) + # Use hollow circles
         geom_point() +
         geom_smooth(method = lm) +
         ggtitle(paste0("Phenotype vs genotype similarity: ", colnames(V), ", heri:", heritability))
    
    
    ggsave(filename = paste0("results/plots/", colnames(V), "_", datetime.stamp, ".pdf"), 
           plot = p, width = 17, height = 7)
    
    # Return
    list(heritability = 0)
    
}

