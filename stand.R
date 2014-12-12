

#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
chr<-para$chr[index]

 testped = read.table(paste("~/NFG/raw/Hd",chr,"clA.raw",sep=""), header=T,sep='')
##testpheno= read.table('~/data/raw/Pheno', header=TRUE, sep='',skip=11)
## gsub("X","NA",testpheno)
## write.table('~/data/raw/pheno2')
 testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
 testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
library(plyr)
 testfn<-testfi
 for (j in 2:ncol(testfi)){is.na(testfn[,j]) <- testfn[,j]=="X"}
 testfnr<-testfn[,c(-15,-17,-18,-31)]
 for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j])
                               testfnr[!ok,j]<- mean(testfnr[ok,j])  }
 testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
 testfrm<-na.omit(testfrm)



snpname= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)

write.table(testfrm,file=paste("~/NFG/raw/leHd",chr,"cl",sep=""))
testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))
## gfrm<-grep("cnv",names(testfrm))
## testfrm<-testfrm[,-gfrm]

dnul<-rep(1==0,ncol(testfrm))
for (j in 13:ncol(testfrm))
{ dnul[j]<-max(testfrm[ ,j])==0}
bnul<-dnul[-(1:13)]
testfrm<-testfrm[,!dnul]
snpnamesous<-snpnamesous[!bnul,]

##compare size in case
ncol(testfrm)
ncol(snpnamesous)

ptdv1=rep(0, (ncol(testfrm)-13) )
for (j in 1:(ncol(testfrm)-13) )
{
if(max(testfrm[,j+13])==0) {ptdv1[j]<-NA} else 
{ptdv1[j]<-summary(lm(testfrm[,4]~testfrm[,j+13]))$coef[2,4]}}


ptdv2=rep(0,(ncol(testfrm)-13) )
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



snpnamesous<-as.data.frame(snpnamesous)
ptdv<-cbind(snpnamesous,ptdv1,ptdv2,ptdv3,ptdv4)

snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")

colnames(ptdv)<-c("snp","KolH","KolL","Trigl","BMI")
ptdv<-join(snpposi[,c(2,4)],ptdv, "snp",type='right')
write.table(ptdv,file=paste("~/NFG/result/Hdapstand",chr,"pv",sep=""))


for (pheno in 3:6 )
{
pdf(file=paste("~/NFG/result/Hstand",names(ptdv)[pheno],"-",chr,"pv",".pdf",sep=""))

plot(ptdv[,2]/1000000,-log(ptdv[,pheno],10) ,col="red",pch=18,xlim=c(40,60),ylim=c(0,13),ann=FALSE)
# points(-log(ptdv[,3],10),col="blue",pch=19)
# points(-log(ptdv[,5],10),col="orange",pch=17)
# points(-log(ptdv[,6],10),col="pink",pch=16)
title(xlab=paste("Chr",chr, "position(Mb) by stand","pheno=",pheno,sep""), col.lab=rgb(0,0.5,0))
title(ylab=paste(names(ptdv)[pheno],":-log10(p-vlue)",sep=""),col.lab=rgb(.5,0,0))


# legend(1,5,c("H","L","tri"),cex=0.8,col=c("blue","orange","pink"),pch=c(19,17,16), lty=1:2)
dev.off()
}