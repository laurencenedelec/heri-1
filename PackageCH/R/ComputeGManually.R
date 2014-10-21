library(doMC)
library(foreach)
registerDoMC(16)

N <- 20

# Get results (stored into SNPS_datetime.stamp.raw)
Phenotypes <- read.table("data/SNPs_14102014_145630.raw", header = TRUE)
# Remove rows with NA's 
Phenotypes <- Phenotypes[complete.cases(Phenotypes), ]
# We must have "famid" and "id" as the two first columns
colnames(Phenotypes) <- c("famid", "id")

# Value of SNPs are stored in the last N columns
snps <- Phenotypes[, (ncol(Phenotypes)-N+1):ncol(Phenotypes)]

build_matrix_G <- function(snps) {
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

#G <- build_matrix_G(snps)
#save(list = "G", file = "G.RData")

