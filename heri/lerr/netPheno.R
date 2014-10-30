
testpheno= read.table('~/NFG/raw/NFBC/NFBC/Subject.Phenotype.txt', header=TRUE, sep='',skip=20)
##  nettoyer le nom 
testpheno[1,]
gsub("X","NA",testpheno)
write.table('~/data/raw/Pheno3')
ncol(testpheno)