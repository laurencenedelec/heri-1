testped = read.table('~/data/test6.raw', header=T,sep='')
testpheno= read.table('~/data/pheno2', header=T, sep='')
testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
testfn<-testfi
for (j in 2:ncol(testfi)){
			is.na(testfn[,j]) <- testfn[,j]=="X"}
			testfnr<-testfn[,c(-15,-17,-18,-31)]
for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
colSums(is.na(testfnr)) 
testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
testfrm<-na.omit(testfnr)
nrow(testfrm)
d<-15
f<-130
pas<-21
dvc=matrix(rep(0,pas),ncol=pas)
dc=matrix(rep(0,pas*9),ncol=9)
tdv=matrix(rep(0,pas*9),ncol=9)
i<-0 
for (k in c(4,5,6,9))
 {d<-14
 i<-i+1
 f<-130
 for (j in 1:pas){ 
     	   	   dc[j,i]<-dcor(testfrm[,k],testfrm[,d:f],1)
		   tdv[j,i]<-dcor.t(testfrm[,k],testfrm[,d:f])
         	   d<-d+100
         	   f<-f+100}
 }
d<-14
 for (j in 1:pas){ dvc[j]<-dcor(testfrm[,d],testfrm[,d+1],1)
           	   d<-d+100}
write.table(dc,file="~/data/disCor100")
d<-dcor(testfrm[,30],testfrm[,31:33])
write.table(dvc,file="~/data/dcorvar")
write.table(tdv,file="~/data/dcort100")

