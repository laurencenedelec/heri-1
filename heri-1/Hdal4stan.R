# testped = read.table('~/NFG/raw/Hdap4.raw', header=T,sep='')
# testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
# testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
library(plyr)
# testfn<-testfi
# for (j in 2:ncol(testfi)){
#			is.na(testfn[,j]) <- testfn[,j]=="X"}
# testfnr<-testfn[,c(-15,-17,-18,-31)]
# for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
#     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
# colSums(is.na(testfnr)) 
# testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
# testfrm<-na.omit(testfrm)
# write.table(testfrm,file="~/data/leHdap4")
testfrm=read.table('~/NFG/raw/leHdap4',sep='')
snpname=read.table('~/NFG/raw/Hdap4.snplist',sep='')
souslist=grep("cnv",names(testfrm))
testfrm<-testfrm[,-souslist]
soussnp=grep("cnv",snpname[,1])
snpnamesous<-snpname[-soussnp,]

<<<<<<< HEAD
# 4=H 5=L 6=T
=======

>>>>>>> d0baed1871cc41ba96847d223cfc63499a75bec3

ptdv1=rep(0, (ncol(testfrm)-13) )
for (j in 1: (ncol(testfrm)-13) )
{var<-testfrm[,j+13]
if(max(var)==0) {ptdv1[j]<-NA} else 
{ptdv1[j]<-summary(lm(testfrm[,4]~var))$coef[2,4]}}


ptdv2=rep(0,ncol(testfrm)-13)
for (j in 1:(ncol(testfrm)-13) )
{var<-testfrm[,j+13]
if(max(var)==0) {ptdv2[j]<-NA} else 
{ptdv2[j]<-summary(lm(testfrm[,5]~var))$coefficients[2,4]}}

ptdv3=rep(0,ncol(testfrm)-13)
for (j in 1:(ncol(testfrm)-13))
{var<-testfrm[,j+13]
if(max(var)==0) {ptdv3[j]<-NA} else {
 ptdv3[j]<-summary(lm(testfrm[,6]~var))$coefficients[2,4]} }

ptdv4=rep(0,ncol(testfrm)-13)
for (j in 1:(ncol(testfrm)-13))
{var<-testfrm[,j+13]
if(max(var)==0) {ptdv4[j]<-NA} 
else {ptdv4[j]<-summary(lm(testfrm[,9]~var))$coefficients[2,4]} }




# for (k in c(4,5,6,9))
# {d<-14
# ki<-ki+1
# f<-14+pas
# for (j in 1:npas)
#     {
 #    dcto[j,ki]<-dcor(testfrm[,k],testfrm[,d:fin],1)    	
  #   tempcor<-dcor.ttest(testfrm[,k],testfrm[,d:f],distance=FALSE)
  #   dc[j,ki]<-tempcor$estimate
  #   ptdv[j,ki]<-tempcor$p.value
  #   Ttdv[j,ki]<-tempcor$statistic
  #   d<-d+pas
  #   f<-f+pas}
 #  }


snpposi<-read.table('~/NFG/raw/snpposi.bim',sep='')
# snpname<-read.table('~/NFG/raw/Hdap4.snplist',sep='')
snpnamesous<-as.data.frame(snpnamesous)
ptdv<-cbind(snpnamesous,ptdv1,ptdv2,ptdv3,ptdv4)
colnames(ptdv)<-c("snp","KolH","KolL","trigl","BMI")
names(snpposi)<-c("x","snp","y","posi","w","z")
ptdv<-join(snpposi[,c(2,4)],ptdv, "snp",type="right")
names(ptdv)<-c("snp","posi","KolH","KolL","trigl","BMI")
write.table(ptdv,file="~/NFG/result/Hdap4pstand")



