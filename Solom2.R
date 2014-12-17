#testped = read.table('~/data/test6.raw', header=T,sep='')
#testpheno= read.table('~/data/pheno2', header=T, sep='')
#testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
#testfn<-testfi
#for (j in 2:ncol(testfi)){
#			is.na(testfn[,j]) <- testfn[,j]=="X"}
#testfnr<-testfn[,c(-15,-17,-18,-31)]
#for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
#     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
#colSums(is.na(testfnr)) 
#testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
#testfrm<-na.omit(testfrm)
#write.table(testfrm,file="~/data/letest")
testfrm=read.table('~/data/letest',sep='')
testfrm[1,1:30]
d<-14
pas<-30
f<-d+pas
npas<-floor((ncol(testfrm)-14)/pas)-2
npas4<-npas*4
fin<-ncol(testfrm)
dvc=matrix(rep(0,npas),ncol=1)
dc=matrix(rep(0,npas4),ncol=4)
tdv=matrix(rep(0,npas4),ncol=4)
ptdv=matrix(rep(0,npas4),ncol=4)
Ttdv=matrix(rep(0,npas4),ncol=4)
dcto=matrix(rep(0,npas4),ncol=4)
i<-0 
#d<-dcor.ttest(testfrm[,30],testfrm[,31:33])
#Ttdv[1,1]<-d$statistic
#Ttdv[1,1]
ki<-0
for (k in c(4,5,6,9))
{d<-14
ki<-ki+1
f<-14+pas
for (j in 1:npas)
     {
     dcto[j,ki]<-dcor(testfrm[,k],testfrm[,d:fin],1)    	
     tempcor<-dcor.ttest(testfrm[,k],testfrm[,d:f],distance=FALSE)
     dc[j,ki]<-tempcor$estimate
     ptdv[j,ki]<-tempcor$p.value
     Ttdv[j,ki]<-tempcor$statistic
     d<-d+pas
     f<-f+pas}
 }
write.table(ptdv,file="~/data/allpvalue30")
write.table(Ttdv,file="~/data/allTvalue30")
write.table(dc,file="~/data/alldc30")
write.table(dcto,file="~/data/alldcto30")
jpeg(file="~/data/allplotp30.jpg")
plot(ptdv[,2])
dev.off()
#d<-14
#for (j in 1:npas){ dvc[j]<-dcor(testfrm[,d],testfrm[,d+1],1)
#    d<-d+pas}
#write.table(dvc,file="~/data/dcorvar10")



