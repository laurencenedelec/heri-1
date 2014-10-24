library(testthat)
##########################
### Load Phenotypes
##########################
P.raw <- read.csv(file = "data/PhenoData.csv")
#load("data/Phenotype_Simulated.RData")
#load("data/Phenotype_Simulated_14102014_104636.RData")
#load("data/Phenotype_Simulated_14102014_145630.RData")
load("data/Phenotype_Simulated_24102014_104043.RData")

P.raw$full_id <- paste(P.raw$famid, P.raw$id)
rownames(P.raw) <- P.raw$full_id
#colnames(P.raw) <- sub("resIN_", "", colnames(P.raw))

# A few tests to assess that our data are correct
#expect_that(dim(P.raw), equals(c(1430, 176)))
#expect_that(colnames(P.raw)[1:3], equals(c("famid", "id", "Affected")))
