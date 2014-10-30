

j<-
pas<-
chr<-19
Repli<-10^10
pheno<-5
library(energy)
print("p permutation")
snpname= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)

testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))

d<-dcov.test(testfrm[,pheno],testfrm[,j:(j+pas)],R=Repli)$p.value
print(d)
print(chr)
print("LDL")
