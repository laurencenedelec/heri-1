# testped = read.table('~/NFG/raw/all19.raw', header=T,sep='')
# testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
# testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(energy)
library(plyr)
# testfn<-testfi
# for (j in 2:ncol(testfi)){
#			is.na(testfn[,j]) <- testfn[,j]=="X"}
#  testfnr<-testfn[,c(-15,-17,-18,-31)]
#  for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
#     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
# colSums(is.na(testfnr)) 
# testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
#  testfrm<-na.omit(testfrm)
# write.table(testfrm,file="~/NFG/raw/leall19")
 testfrm=read.table('~/NFG/raw/leall19',sep='')
snpname=read.table('~/NFG/raw/all19.snplist',sep='')
souslist=grep("cnv",names(testfrm))
testfrm<-testfrm[,-souslist]
soussnp=grep("cnv",snpname[,1])
snpnamesous<-snpname[-soussnp,]



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





snpposi<-read.table('~/NFG/raw/snpposi.bim',sep='')
snpnamesous<-as.data.frame(snpnamesous)
ptdv<-cbind(snpnamesous,ptdv1,ptdv2,ptdv3,ptdv4)
colnames(ptdv)<-c("snp","KolH","KolL","trigl","BMI")
names(snpposi)<-c("x","snp","y","posi","w","z")
ptdv<-join(snpposi[,c(2,4)],ptdv, "snp",type="right")
names(ptdv)<-c("snp","posi","KolH","KolL","trigl","BMI")
write.table(ptdv,file="~/NFG/result/full19pstand")



