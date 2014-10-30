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
    
    # Remove diagonal of Y.square and phi
    diag(Y.square) <- NA
    diag(phi) <- NA
    
    # Recast matrix to have the correct form
    Y.square <- t(matrix(Y.square[which(!is.na(Y.square))],
                         nrow = (nrow(Y.square)-1),
                         ncol = ncol(Y.square)))
    phi <- t(matrix(phi[which(!is.na(phi))],
                    nrow = (nrow(phi)-1),
                    ncol = ncol(phi)))    
    
    # We need this as a vector
    Y.square.vec <- as.vector(Y.square)
    phi <- as.vector(phi)
    
    # We are creating three plots
    # 1) all data
    # 2) only "independent individuals" (< 0.15 genotype similarity)
    # 3) only "dependent individuals" (> 0.15 genotype similarity)
    create_ggplot <- function(df, V, percentage_points_shown = 0.05) {
        ggplot(data = df[sample(1:nrow(df), percentage_points_shown * nrow(df), replace=FALSE),],
                    aes(y = phenotype.similarity, 
                        x = genotype.similarity)) + # Use hollow circles
            geom_point() +
            geom_smooth(method = lm) +
            ggtitle(paste0("Phenotype vs genotype similarity: ", colnames(V), ", 
                        heri:", summary(lm(df$phenotype.similarity ~ df$genotype.similarity))$coef[2,1]))
    }
    
    df <- data.frame(phenotype.similarity = Y.square.vec,
                     genotype.similarity = phi)
    p <- create_ggplot(df, V)

    threshold <- 0.15
    df <- data.frame(phenotype.similarity = Y.square.vec[which(phi >= threshold)],
                     genotype.similarity = phi[which(phi >= threshold)])
    q <- create_ggplot(df, V, percentage_points_shown = 0.5)
    
    df <- data.frame(phenotype.similarity = Y.square.vec[which(phi < threshold)],
                     genotype.similarity = phi[which(phi < threshold)])
    r <- create_ggplot(df, V, percentage_points_shown = 0.025)
    
    # Plot results, 3 columns (= 3 methods)
    pdf(file = paste0("results/plots/filtered_", colnames(V), "_", datetime.stamp, ".pdf"), 
        width = 17, 
        height = 7)
    
        # Create a grid with 2 rows, length(p) columns
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(1, 3)))   

            print(p, vp = vplayout(1,1))
            print(q, vp = vplayout(1,2))        
            print(r, vp = vplayout(1,3))
    
        upViewport(0)
    
    dev.off()
    
    # Return
    list(heritability = 0)
    
}

