library(doMC)
library(foreach)
registerDoMC(16)

setwd("/home/duvanel/git/PackageCH/")

build_matrix_G <- function(path.snps = "data/SNPs_14102014_145630.raw") {
    
    N <- 20
    
    # Get results (stored into SNPS_datetime.stamp.raw)
    Phenotypes <- read.table(path.snps, header = TRUE)
    # Remove rows with NA's 
    Phenotypes <- Phenotypes[complete.cases(Phenotypes), ]
    # We must have "famid" and "id" as the two first columns
    colnames(Phenotypes) <- c("famid", "id")
    
    # Value of SNPs are stored in the last N columns
    snps <- Phenotypes[, (ncol(Phenotypes)-N+1):ncol(Phenotypes)]
    
    G <- matrix(0, 
                ncol = nrow(snps),
                nrow = nrow(snps))
    
    r <- foreach(i=1:nrow(snps)) %dopar% {
        for(j in 1:nrow(snps)) {
            
            if(j %% 10 == 0) cat("i = ", i, "/", nrow(snps), " and j = ", j, "/", nrow(snps), "\n")
            
            G[i,j] <- sqrt(sum((snps[i,] - snps[j,]) ^ 2))
        }
    }
    G
} 

HeritabilityEstimation <- function() {
    
    G <- build_matrix_G(snps)
    load("data/Phenotype_Simulated_14102014_145630.RData")
    
    P <- as.matrix(P.raw$fake.trait)
    
    ### We keep one phenotype and normalize its
    V <- matrix(P, 
                nrow = nrow(P),
                ncol = ncol(P))
    
    V <- t(t(V) - colMeans(V))
    V.sd <- colSds(V)
    
    # Return
    V <- t(t(V) * (1/V.sd))
    
    # Get data for the specific phenotype   
    Y.distance.carre <- matrix(rep(V^2,nrow(V)),
                               ncol = nrow(V)) + 
        t(matrix(rep(V^2,nrow(V)),
                 ncol = nrow(V))) - 
        2 * V %*%t(V)
    
    Y.distance <- sqrt(Y.distance.carre)
    
    heritability <- dcov(Y.distance, G) / (sqrt(dcov(G, G)) * sqrt(dcov(Y.distance, Y.distance)))
    
}
