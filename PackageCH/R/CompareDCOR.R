library(energy)
library(matrixStats)
library(reshape2)

#' Build alpha s.t. X = G * alpha
#'
#' @title Build a random SNps matrix
#' @param N_real_coeff number of non_zero coeff
#' @param N_SNPs number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha
#' @author Julien Duvanel
#' @export
build_alpha <- function(N_real_coeff, N_SNPS, b) {
    
    # build alpha
    alpha <- sample(c(runif(N_real_coeff, b[1], b[2]), 
                      rep(x = 0, times = N_SNPS - N_real_coeff)))
    
}

#' Build a random SNps matrix
#'
#' @title Build a random SNps matrix
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
#' @return export a pdf file
#' @author Julien Duvanel
#' @export
compare_dcor <- function(N, 
                         N_SNPS, 
                         N_real_coeff, 
                         get_snps_matrix,
                         get_alpha,
                         b,
                         variable = "") {
    
    # We have to check that dimensions agree
    if(length(N) != length(N_SNPS) | length(N) != length(N_real_coeff)) stop("Problem, length of N has to be the same as N_SNPS.")    
    if(length(b) != 2) stop("b has to be a vector of length 2 because we use it as runif(..., b[1], b[2]) !")
    
    # We loop through length(N) (but they all have same length at this point)
    res <- c()
    for(i in 1:length(N)) {
        
        # Build a fake genome matrix
        M <- get_snps_matrix(N[i], N_SNPS[i])
        alpha <- get_alpha(N_real_coeff = N_real_coeff[i], 
                           N_SNPS = N_SNPS[i], 
                           b = b)
                
        # Build fake trait X
        X <- M %*% alpha

        # Do dcor estimation
        res <- rbind(res, 
                     c(get(variable[1])[i], 
                       dcor(X, M),
                       dcor(X + rnorm(N[i], mean = 0, sd = 1), M),
                       # This last line is used to compare dcor(X,Y)
                       # when X and Y are completly independent
                       # (this is done by simulating X and then generating a new SNPs matrix)
                       dcor(X, get_snps_matrix(N[i], N_SNPS[i]))))
        
        cat("-> i = ", i, "/", length(N), "\n")
    }
    
    res <- data.frame(var = res[,1],
                      X_from_G = res[,2],
                      X_from_G_plus_noise = res[,3],
                      X_not_from_G = res[,4])

    # Export a pdf file
    p <- ggplot(data = melt(res, measure.vars = c("X_from_G", "X_from_G_plus_noise", "X_not_from_G")),
                aes_string(x = "var" , y = "value")) +
            geom_point(aes_string(color = "variable"), size = 2, position = position_jitter(w = 1.5, h = 0)) +
            xlab(paste0("Value of ", paste(variable, collapse=", "))) +
            ylab("dcor(X,G)") + 
            ggtitle(paste0("X = G * vec of runif(", b[1], ",", b[2],  "), N in ", min(N), ":", max(N),
                           ", N_SNPS in ", min(N_SNPS), ":", max(N_SNPS),
                           " and N_real_coeff in ", min(N_real_coeff), ":", max(N_real_coeff)))

    ggsave(plot = p,
           filename = paste0("results/plots/ExperimentalDcor_", 
                             format(Sys.time(), "%d%m%Y_%H%M%S"),
                             ".pdf"),
           width = 17, height = 7)
    
    save(list = "res", file = paste0("results/ExperimentalDcor_", 
                                     format(Sys.time(), "%d%m%Y_%H%M%S"),
                                     ".RData"))
    # return
    list(res = res, X = X, M = M)    
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
