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
####  load Data for Theoretical Kinship
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
load_data_TheoKin <- function() {

  load("data/TheoKin/K.matrix.RData")
  
  # reorder column
  phi.matrix <- phi.matrix[sort(rownames(phi.matrix)), sort(rownames(phi.matrix))]  
  
  list(phi = phi.matrix)
                          
}

build_K_TheoKin_Once <- function() {
  
  K.raw <- read.csv(file = "data/TheoKin/TheoKinship.txt", sep = " ")
  
  # Create a column with full ID (which is unique)
  K.raw$full_ID1 <- paste(K.raw$Fam, K.raw$PID1)
  K.raw$full_ID2 <- paste(K.raw$Fam, K.raw$PID2)
  colnames(K.raw)[4] <- "value"
  
  # A few tests to assess that our data are correct
  expect_that(dim(K.raw), equals(c(60313, 6)))
  expect_that(colnames(K.raw), equals(c("Fam", "PID1", "PID2", "value", "full_ID1", "full_ID2")))
  
  # We build the matrix phi filling
  # every single connection we have
  # (later we'll filter them)
  phi.matrix <- as.matrix(2 * build_matrix_K(K.raw)$K.matrix)
  
  save(list = "phi.matrix", file="data/TheoKin/K.matrix.RData")
}
