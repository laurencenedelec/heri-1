testped = read.table('~/data/test6.raw', header=T,sep='')
testpheno= read.table('~/data/pheno2', header=T, sep='')
testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
testfn<-testfi
for (j in 2:ncol(testfi)){is.na(testfn[,j]) <- testfn[,j]=="X"}
testfnr<-testfn[,c(-15,-17,-18,-31)]
testfrm<-na.omit(testfnr)
#d<-30
#f<-130
#pas<-22
#dc<-cbind(rep(0,pas),rep(0,pas),rep(0,pas),rep(0,pas),rep(0,pas),rep(0,pas),rep#(0,pas),rep(0,pas),rep(0,pas))
#i<-0 
#for (k in c(5,7,8,10,11,12,22,23,24))
# {d<-30 
# i<-i+1
# f<-130
#     for (j in 1:pas){ dc[j,i]<-dcor(testfrm[,k],testfrm[,d:2270],1)
#         d<-d+100
#         f<-f+100}
# }
#dc

#save.image() 
#sink()
u<-0
dtoto<-dcor.ttest(testfrm[,33],testfrm[,34],distance=FALSE)
dtoto$p.value

v<-dtoto$estimate
v
