library(energy)
library(matrixStats)
library(reshape2)

#' Build alpha s.t. X = G * alpha
#'
#' @title Build alpha
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

#' Build alpha s.t. X = G * alpha but alpha(g_i)
#'
#' @title Build alpha
#' @param N_real_coeff number of non_zero coeff
#' @param N_SNPs number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha with dominant effect
#' @author Julien Duvanel
#' @export
build_alpha_dominant <- function(N_real_coeff, N_SNPS, b) {
    
    # We create an alpha for the dominant effect
    alpha_d <- matrix(0, nrow = 3, ncol = N_SNPS)
    alpha_d[1, ] <- sample(c(runif(N_real_coeff, -1, 1), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    alpha_d[2, ] <- sample(c(runif(N_real_coeff, 3, 4), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    alpha_d[3, ] <- sample(c(runif(N_real_coeff, 15, 20), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    
    alpha_d

}

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for dominant effect
#' @param M a matrix with N_SNPS columns
#' @param alpha a matrix with 3 rows and N_SNPS columsn
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha_dominant <- function(M, alpha) {
    # Usually, as.matrix(mat) %*% alpha does the job
    # but here, the multiplication differs if mat[i,j] has 
    # a specific value.
    
    # For every row, we do the "special product"
    res <- apply(M, 1, function(x) {
        ret <- 0
        for(i in 1:length(x)) {
            ret <- ret + x[i]*alpha[x[i]+1, i]
        }
        ret
    })
    
    as.matrix(res)
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
                         variable = "",
                         product_snps_alpha = function(M, alpha) { M %*% alpha }) {
    
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
        X <- product_snps_alpha(M, alpha)

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
            scale_y_continuous(limits = c(0, 1)) +
            ggtitle(paste0("X = G * vec of runif(", b[1], ",", b[2],  "), N in ", min(N), ":", max(N),
                           ", N_SNPS in ", min(N_SNPS), ":", max(N_SNPS),
                           " and N_real_coeff in ", min(N_real_coeff), ":", max(N_real_coeff)))

    ggsave(plot = p,
           filename = paste0("results/plots/ExperimentalDcor_", 
                             as.character(substitute(get_alpha)), "_",
                             format(Sys.time(), "%d%m%Y_%H%M%S"),
                             ".pdf"),
           width = 17, height = 7)
    
    save(list = "res", file = paste0("results/ExperimentalDcor_", 
                                     as.character(substitute(get_alpha)), "_",
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
