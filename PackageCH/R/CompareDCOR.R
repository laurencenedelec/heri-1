library(energy)
library(matrixStats)

#' Build a random SNps matrix
#'
#' @title GBuild a random SNps matrix
#' @param N sample size
#' @param N_SNPs number of SNPS
#' @param snps_value value that SNPs can take
#' @return a random SNPs matrix
#' @author Julien Duvanel
#' @export
build_SNPs_matrix <- function(N, N_SNPS, snps_value = c(0,1,2)) {

    # build a completly random matrix
    M <- matrix(sample(snps_value, size = N * N_SNPS, replace = T), 
                nrow = N)
    
}

#' Compare DCOR with different kind of settings
#'
#' @title GBuild a random SNps matrix
#' @param N sample size
#' @param N_SNPs number of SNPS
#' @param N_real_coeff number of SNPs that really explain the phenotype
#' @param get_snps_matrix function that gives a SNPs matrix
#' @param b X = G * runif(, -b, b)
#' @param variable the variable that vary
#' @param AddNoise if T it adds noise to X = G * alpha + epsilon with N(0,1)
#' @return export a pdf file
#' @author Julien Duvanel
#' @export
compare_dcor <- function(N, 
                         N_SNPS, 
                         N_real_coeff, 
                         get_snps_matrix,
                         b,
                         variable = "",
                         AddNoise = FALSE) {
    
    # We have to check that dimensions agree
    if(length(N) != length(N_SNPS) | length(N) != length(N_real_coeff)) stop("Problem, length of N has to be the same as N_SNPS.")    
    if(length(b) != 2) stop("b has to be a vector of length 2 because we use it as runif(..., b[1], b[2]) !")
    
    # We loop through length(N) (but they all have same length at this point)
    res <- c()
    for(i in 1:length(N)) {
        
        # Build a fake genome matrix
        M <- get_snps_matrix(N[i], N_SNPS[i])
        alpha <- sample(c(runif(N_real_coeff[i], b[1], b[2]), rep(x = 0, times = N_SNPS[i] - N_real_coeff[i])))
        
        # Build fake trait X
        X <- M %*% alpha
        
        if(AddNoise == TRUE) {
            X <- X + rnorm(N[i], mean = 0, sd = 0)
        }
        
        # Do dcor estimation
        res <- rbind(res, 
                     c(get(variable[1])[i], dcor(X, M)))
        
        cat("-> i = ", i, "/", length(N), "\n")
    }
    
    # Export a pdf file
    pdf(file = paste0("results/plots/ExperimentalDcor_", format(Sys.time(), "%d%m%Y_%H%M%S") ,".pdf"), width = 17, height = 7)
    
        txtNoise <- ""
        if(AddNoise == TRUE) txtNoise <- " + N(0,1)"
        
        plot(res, 
             xlab = paste0("Value of ", paste(variable, collapse=", ")), 
             ylab = "dcor(X,G)",
             main = paste0("X = G * vec of runif(", b[1], ",", b[2],  ")", txtNoise, ", N in ", min(N), ":", max(N),
                           ", N_SNPS in ", min(N_SNPS), ":", max(N_SNPS),
                           " and N_real_coeff in ", min(N_real_coeff), ":", max(N_real_coeff)))
    
    dev.off()
    
    list( res = res)
    
}






# X.cov <- cov(X)
# M.cov <- cov(M)
# 
# X.tilde <- X %*% sqrtm(solve(X.cov))
# M.tilde <- M %*% sqrtm(solve(M.cov))
# dcor(X.tilde, M.tilde)



# 
# N <- 200
# res <- c()
# for(i in 1:200) {
#     
#     a <- c(0,1,2)
#     N <- 20
#     N_SNPS <- 20
#     M <- matrix(sample(a, size = N * N_SNPS, replace = T), 
#                 nrow = N)
#     alpha <- runif(N_SNPS, 4, 4)
#     
#     X <- M %*% alpha
#     
#     #X <- rnorm(n = N, mean = 0, sd = 1)
#     #Y <- rnorm(n = N, mean = 0, sd = 1)
#     
#     res <- rbind(res, cbind(
#                             dcor(X, M)^2,
#                             
#                             abs(cov(X,M)/var(M)))
#                             )
# }
# plot(res)
# lm.res <- lm(res[,2] ~ res[,1])
# 
