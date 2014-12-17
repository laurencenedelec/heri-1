#testped = read.table('~/data/raw/Hap4.raw', header=T,sep='')
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
testfrm=read.table('~/data/raw/letest',sep='')
testfrm[1,1:30]
d<-14
pas<-30
f<-d+pas
npas<-floor((ncol(testfrm)-14)/pas)-pas-1
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

d<-14
f<-14+pas
for (j in 1:npas)
     {
     #dcto[j]<-dcor(testfrm[,c(4,5,6,7)],testfrm[,d:fin],1)    	
     tempcor<-dcor.ttest(testfrm[,c(4,5,6,7)],testfrm[,d:f],distance=FALSE)
     dc[j]<-tempcor$estimate
     ptdv[j]<-tempcor$p.value
     Ttdv[j]<-tempcor$statistic
     d<-d+pas
     f<-f+pas}
 
write.table(ptdv,file="~/data/Mallpvalue30")
write.table(Ttdv,file="~/data/MallTvalue30")
write.table(dc,file="~/data/Malldc30")
#write.table(dcto,file="~/data/Malldcto30")
jpeg(file="~/data/Mallplotp30.jpg")
plot(ptdv[])
dev.off()




