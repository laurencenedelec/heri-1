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
####  load Data for GCTA
####
####**********************************************************************

# library used in this R script
library()

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
load_data_IBD <- function() {

  load("data/IBD/phi.matrix.RData")
  
  list(phi = phi.matrix)
                          
}

build_phi_matrix_once <- function() {
  phi.matrix <- read.table("data/IBD/plink.genome", header = TRUE)
  phi.matrix$full_ID1 <- paste(phi.matrix$FID1, phi.matrix$IID1)
  phi.matrix$full_ID2 <- paste(phi.matrix$FID2, phi.matrix$IID2)
  
  phi.temp <- phi.matrix[, c("full_ID1", "full_ID2", "EZ")]
  colnames(phi.temp) <- c("full_ID1", "full_ID2", "value")
  
  phi.matrix <- 0.5 * build_matrix_K(K = phi.temp)$K.matrix
}