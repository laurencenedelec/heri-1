library(testthat)
##########################
### Load Phenotypes
##########################
P.raw <- read.csv(file = "data/PhenoData.csv")
P.raw$full_id <- paste(P.raw$famid, P.raw$id)
rownames(P.raw) <- P.raw$full_id
colnames(P.raw) <- sub("resIN_", "", colnames(P.raw))

# A few tests to assess that our data are correct
expect_that(dim(P.raw), equals(c(1430, 176)))
expect_that(colnames(P.raw)[1:3], equals(c("famid", "id", "Affected")))