library(energy)
library(matrixStats)

build_SNPs_matrix <- function(N, N_SNPs, snps_value = c(0,1,2)) {
    
    M <- matrix(sample(snps_value, size = N * N_SNPS, replace = T), 
                nrow = N)
    
}

compare_dcor <- function(N, 
                         N_SNPS, 
                         N_real_coeff, 
                         get_snps_matrix,
                         b,
                         variable = "") {
    
    if(length(N) != length(N_SNPS) | length(N) != length(N_real_coeff)) stop("Problem, length of N has to be the same as N_SNPS.")    
    if(length(b) != 2) stop("b has to be a vector of length 2 because we use it as runif(..., b[1], b[2]) !")
    
    res <- c()
    for(i in 1:length(N)) {
        
        # Build a fake genome matrix
        M <- get_snps_matrix(N[i], N_SNPS[i])
        alpha <- sample(c(runif(N_real_coeff[i], b[1], b[2]), rep(x = 0, times = N_SNPS[i] - N_real_coeff[i])))
        
        # Build fake trait X
        X <- M %*% alpha
        
        res <- rbind(res, 
                     c(get(variable)[i], dcor(X, M)))
        
        cat("i = ", i, "\n")
    }
    
    pdf(file = paste0("results/plots/ExperimentalDcor_", format(Sys.time(), "%d%m%Y_%H%M%S") ,".pdf"), width = 17, height = 7)
    
    plot(res, 
         xlab = paste0("Value of ", variable), 
         ylab = "dcor(X,G)",
         main = paste0("X = G * vec of runif(", b[1], ",", b[2],  "), N in ", min(N), ":", max(N),
                       ", N_SNPS in ", min(N_SNPS), ":", max(N_SNPS),
                       " and N_real_coeff in ", min(N_real_coeff), ":", max(N_real_coeff)))
    
    dev.off()
    
}

N <- c(100, 1000)
N_SNPS <- c(100, 100)
N_real_coeff <- c(20, 20)
b <- c(2,2)

compare_dcor(N, N_SNPS,N_real_coeff, build_SNPs_matrix, b, "N")






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
