# testped = read.table('~/data/test6.raw', header=T,sep='')
# testpheno= read.table('~/data/pheno2', header=T, sep='')
# testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
# testfn<-testfi
# for (j in 2:ncol(testfi)){
#			is.na(testfn[,j]) <- testfn[,j]=="X"}
# testfnr<-testfn[,c(-15,-17,-18,-31)]
# for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
#     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
# colSums(is.na(testfnr)) 
# testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
# testfrm<-na.omit(testfrm)
# write.table(testfrm,file="~/data/letest")
testfrm=read.table('~/NFG/raw/leHdap19',sep='')
#testfrm[1,1:30]

dc=matrix(rep(0,3000*100),ncol=100)
## ptdv=matrix(rep(0,3000*100),ncol=100)
## Ttdv=matrix(rep(0,3000*100),ncol=100)

for (pas in 1:50) 
{d<-14
k<-5
f<-d+pas
npas<-pmax(floor((ncol(testfrm)-14)/pas)-2,1)
for (j in 1:npas)
     {
    ## tempcor<-dcor.ttest(testfrm[,k],testfrm[,d:f],distance=FALSE)
     dc[j,pas]<-dcor(testfrm[,k],testfrm[,d:f],1)
    ## ptdv[j,pas]<-tempcor$p.value
    ## Ttdv[j,pas]<-tempcor$statistic
     d<-d+pas
     f<-f+pas}
 }
## write.table(ptdv,file="~/data/p5pas")
## write.table(Ttdv,file="~/data/T5pas")
write.table(dc,file="~/NFG/result/tailledcpas")





