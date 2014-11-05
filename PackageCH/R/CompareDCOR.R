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

#' Build alpha s.t. X = G * alpha but alpha(g_i,g_j)
#'
#' @title Build alpha
#' @param N_real_coeff number of non_zero coeff
#' @param N_SNPs number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha with dominant effect
#' @author Julien Duvanel
#' @export 
build_alpha_epistatic <- function(N_real_coeff, N_SNPS, b) {
    
    # We create an alpha for the dominant effect
    alpha_d <- matrix(0, nrow = 4, ncol = N_SNPS)
    
    # If product of g_i and g_j = 0
    alpha_d[1, ] <- sample(c(runif(N_real_coeff, -1, 1), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    # If product of g_i and g_j = 1
    alpha_d[2, ] <- sample(c(runif(N_real_coeff, 3, 4), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    # If product of g_i and g_j = 2
    alpha_d[3, ] <- sample(c(runif(N_real_coeff, 13, 15), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))
    # If product of g_i and g_j = 4
    alpha_d[4, ] <- sample(c(runif(N_real_coeff, 15, 20), 
                             rep(x = 0, times = N_SNPS - N_real_coeff)))    
    alpha_d
}

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for additive effect
#' @param M a matrix with N_SNPS columns
#' @param alpha a vector 
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha <- function(M, alpha) { M %*% alpha }

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

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for epistatic effect
#' @param M a matrix with N_SNPS columns
#' @param alpha a matrix with 3 rows and N_SNPS columsn
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha_epistatic <- function(M, alpha) {
    # Usually, as.matrix(mat) %*% alpha does the job
    # but here, the multiplication differs if mat[i,j] has 
    # a specific value.
    
    f <- function(x) {
        if(x == 0) 1
        else if(x == 1) 2
        else if(x == 2) 3
        else if(x == 4) 4
        else 1
    }
    # For every row, we do the "special product"
    res <- apply(M, 1, function(x) {
        ret <- 0
        for(i in 1:(length(x)-1)) {
            ret <- ret + x[i]*x[i+1]*alpha[f(x[i]*x[i+1]), i]
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
#' @title Compare DCOR and LM
#' @param N sample size
#' @param N_SNPs number of SNPS
#' @param N_real_coeff number of SNPs that really explain the phenotype
#' @param b X = G * runif(, -b, b)
#' @param delta_add importance of additive effect
#' @param delta_dom importance of dominant effect
#' @param delta_epi importance of epistatic effect
#' @param variable the variable that vary
#' @return export a pdf file
#' @author Julien Duvanel
#' @export
compare_dcor <- function(N, 
                         N_SNPS, 
                         N_real_coeff, 
                         b,
                         delta_add,
                         delta_dom,
                         delta_epi,
                         variable = "") {
    
    # We have to check that dimensions agree
    if(length(N) != length(N_SNPS) | length(N) != length(N_real_coeff)) stop("Problem, length of N has to be the same as N_SNPS.")    
    if(length(b) != 2) stop("b has to be a vector of length 2 because we use it as runif(..., b[1], b[2]) !")
    
    # We loop through length(N) (but they all have same length at this point)
    res <- c()
    for(i in 1:length(N)) {
        
        # Build two fake genome matrices
        # the second one is only used to compare with random results
        M <- build_SNPs_matrix(N[i], N_SNPS[i])
        M_tilde <- build_SNPs_matrix(N[i], N_SNPS[i])
        
        # Prepare all alpha's (additive/dominant/epistatic effect)
        alpha_add <- build_alpha(N_real_coeff = N_real_coeff[i], 
                                 N_SNPS = N_SNPS[i], 
                                 b = b)
        alpha_dom <- build_alpha_dominant(N_real_coeff = N_real_coeff[i], 
                                          N_SNPS = N_SNPS[i], 
                                          b = b)
        
        alpha_epi <- build_alpha_epistatic(N_real_coeff = N_real_coeff[i], 
                                           N_SNPS = N_SNPS[i], 
                                           b = b)
        
        # Build fake trait X
        X <- delta_add[i] * product_snps_alpha(M, alpha_add) + 
             delta_dom[i] * product_snps_alpha_dominant(M, alpha_dom) +
             delta_epi[i] * product_snps_alpha_epistatic(M, alpha_epi)
        
        # We are more interested to have X + noise 
        # noise is 10% of sd(X)
        noise <- rnorm(N[i], mean = 0, sd = 0.1 * ifelse(sd(X) == 0, 1, sd(X)))   
        X_norm_plus_noise <- (X + noise - mean(X + noise)) / sd(X + noise)
        
        # Compute in advance distance matrices for
        # linear regression estimate of heritability
        dist_M <- as.vector(as.matrix(dist(M)))
        dist_M_tilde <- as.vector(as.matrix(dist(M_tilde)))

        # Do estimates
        res <- rbind(res, 
                     c(get(variable[1])[i], 

                       # dcor estimates
                       dcor(X_norm_plus_noise, M),
                       dcor(X_norm_plus_noise, M_tilde),
    
                       # lm estimates
                       lm(as.vector(X_norm_plus_noise  %*% t(X_norm_plus_noise)) ~ dist_M)$coefficients[2],
                       lm(as.vector(X_norm_plus_noise %*% t(X_norm_plus_noise)) ~ dist_M_tilde)$coefficients[2]))
        
        cat("-> i = ", i, "/", length(N), "\n")
    }
    
    # Gather data into a dataframe
    res <- data.frame(var = res[,1],
                      
                      dcor_X_from_G_plus_noise = res[,2],
                      dcor_X_not_from_G = res[,3],
                      
                      lm_from_G_plus_noise = res[,4],
                      lm_not_from_G = res[,5])
    
    # melt data to be able plot group into ggplots
    data.melt <- melt(res, measure.vars = c("dcor_X_from_G_plus_noise", 
                                            "dcor_X_not_from_G",
                                            
                                            "lm_from_G_plus_noise",
                                            "lm_not_from_G"))
        
    # Export a pdf file
    p <- ggplot(data = data.melt,
                aes_string(x = "var" , y = "value")) +
            geom_point(aes_string(color = "variable"), size = 2, position = position_jitter(w = 0.005 * sd(data.melt$value), h = 0)) +
            xlab(paste0("Value of ", paste(variable, collapse=", "))) +
            ylab("dcor(X,G)") + 
            ggtitle(paste0("X = G * vec of runif(", b[1], ",", b[2],  "), N in ", min(N), ":", max(N),
                           ", N_SNPS in ", min(N_SNPS), ":", max(N_SNPS),
                           ", delta add in ", min(delta_add), ":", max(delta_add),
                           ", dom in ", min(delta_dom), ":", max(delta_dom),
                           ", epi in ", min(delta_epi), ":", max(delta_epi),
                           " and N_real_coeff in ", min(N_real_coeff), ":", max(N_real_coeff)))

    ggsave(plot = p,
           filename = paste0("results/plots/dcor_vs_lm_", 
                             as.character(substitute(get_alpha)), "_",
                             format(Sys.time(), "%d%m%Y_%H%M%S"),
                             ".pdf"),
           width = 17, height = 7)
    
    save(list = "res", file = paste0("results/data/dcor_vs_lm_", 
                                     as.character(substitute(get_alpha)), "_",
                                     format(Sys.time(), "%d%m%Y_%H%M%S"),
                                     ".RData"))
    
    # return
    list(res = res, X = X, M = M)    
    
}
