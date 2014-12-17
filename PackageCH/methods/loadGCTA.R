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

#' 
#'
#' @title 
#' @param 
#' @return 
#' @author Julien Duvanel
#' @export
load_data_GCTA <- function() {

  data_GCTA <- read.table("data/GCTA/plink.rel")
  data_GCTA_id <- read.table("data/GCTA/plink.rel.id")
  
  colnames(data_GCTA) <- paste(data_GCTA_id$V1, data_GCTA_id$V2)
  rownames(data_GCTA) <- paste(data_GCTA_id$V1, data_GCTA_id$V2)

  # Data are stored in the data/METHODNAME directory
  #data_GCTA <- ReadGRMBin(prefix = paste0("data/GCTA/BP.GRM"))
  
  # Create the matrix \varphi from the data
  # Data come as a lower triangular matrix
  #varphi <- diag(nrow(data_GCTA$id))
  #diag(varphi) <- data_GCTA$diag
  #varphi[lower.tri(varphi, diag = F)] <- data_GCTA$off
  #varphi <- varphi + t(varphi) - diag(diag(varphi))
  
  # Give the correct col/row names (this information is used later on to filter data)
  #colnames(varphi) <- paste(data_GCTA$id$V1, data_GCTA$id$V2)
  #rownames(varphi) <- paste(data_GCTA$id$V1, data_GCTA$id$V2)
  
  list(phi = as.matrix(data_GCTA))
                          
}
