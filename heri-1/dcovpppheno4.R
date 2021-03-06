#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
## dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
pas <- para$pas[index]
chr<-para$chr[index]
Repli<-1000
library(energy)
library(plyr)

testped = read.table(paste("~/NFG/raw/Hd",chr,"clA.raw",sep=""), header=T,sep='')
##testpheno= read.table('~/data/raw/Pheno', header=TRUE, sep='',skip=11)
## gsub("X","NA",testpheno)
## write.table('~/data/raw/pheno2')
 testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
 testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))
 testfn<-testfi
 for (j in 2:ncol(testfi)){is.na(testfn[,j]) <- testfn[,j]=="X"}
 testfnr<-testfn[,c(-15,-17,-18,-31)]
 for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
## colSums(is.na(testfnr)) 
testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
testfrm<-na.omit(testfrm)

snpname= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)

## cut<-( (nrow(snpnamesous)-9):nrow(snpnamesous))
## snpsoussous<-snpnamesous[-cut,]
## snpsoussous<-as.data.frame(snpsoussous)
 
write.table(testfrm,file=paste("~/NFG/raw/leHd",chr,"cl",sep=""))
 
## testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))
## gfrm<-grep("cnv",names(testfrm))
## testfrm<-testfrm[,-gfrm]
 
##compare size in case
ncol(testfrm)
ncol(snpnamesous)

##cleaning the zero column
dnul<-rep(1==0,ncol(testfrm))
for (j in 13:ncol(testfrm)) 
{ dnul[j]<-max(testfrm[ ,j])==0}
bnul<-dnul[-(1:13)]
testfrm<-testfrm[,!dnul]
snpnamesous<-snpnamesous[!bnul,]
snpnamesous<-as.data.frame(snpnamesous)

##compare size in case
ncol(testfrm)
ncol(snpnamesous)

##compute the pvalue for only L
ndonne<-(ncol(testfrm)-13)*4
dc=matrix(rep(0,ndonne),ncol=4)
ptdv=matrix(rep(0,ndonne),ncol=4)
Ttdv=matrix(rep(0,ndonne),ncol=4)
dcsin=matrix(rep(0,ndonne),ncol=4)
dcpsin=matrix(rep(0,ndonne),ncol=4)
dcb=matrix(rep(0,ndonne),ncol=4)

 
for (j in 14:(ncol(testfrm)-pas) )
     {  	
     tempcor<-dcov.test(testfrm[,5],testfrm[,j:(j+pas)],R=Repli)
     dcb[j-13,2]<-tempcor$estimate
     dc[j-13,2]<-dcor(testfrm[,5],testfrm[,j:(j+pas)],1)
     ptdv[j-13,2]<-tempcor$p.value
     Ttdv[j-13,2]<-tempcor$statistic
     }



snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")


ptdv<-cbind(snpnamesous,ptdv)
Ttdv<-cbind(snpnamesous,Ttdv)
dc<-cbind(snpnamesous,dc)
dcb<-cbind(snpnamesous,dcb)



names(ptdv)<-c("snp","KolH","KolL","Trigl","BMI")
names(Ttdv)<-c("snp","KolH","KolL","Trigl","BMI")
names(dc)<-c("snp","KolH","KolL","Trigl","BMI")
names(dcb)<-c("snp","KolH","KolL","Trigl","BMI")


ptdv<-join(snpposi[,c(2,4)],ptdv, 'snp',type='right')
dc<-join(snpposi[,c(2,4)],dc,'snp', type='right')
Ttdv<-join(snpposi[,c(2,4)],Ttdv, 'snp',type='right')
dcb<-join(snpposi[,c(2,4)],dcb,'snp', type='right')

write.table(ptdv,file=paste("~/NFG/result/Hdap",chr,"pv",pas,sep=""))
write.table(Ttdv,file=paste("~/NFG/result/Hdap",chr,"Tv",pas,sep=""))
write.table(dc,file=paste("~/NFG/result/Hdap",chr,"dc",pas,sep=""))
write.table(dcb,file=paste("~/NFG/result/Hdap",chr,"dcb",pas,sep=""))

## remove zero 
## ptdvz<-ptdv[,4]
## ptdvz[ptdvz==0]<-NA

## draw pvalues
pheno<-4  

pdf(file=paste("~/NFG/result/Hd",names(ptdv)[pheno],"-",chr,"pvR",pas,".pdf",sep=""))

plot(ptdv[,2]/1000000,-log(ptdv[,pheno],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,13),ann=FALSE)
# points(-log(ptdv[,3],10),col="blue",pch=19)
# points(-log(ptdv[,5],10),col="orange",pch=17)
# points(-log(ptdv[,6],10),col="pink",pch=16)
title(xlab=paste("Chr",chr, "position(Mb) by dcov-pas=",pas,"pheno=",pheno,,"replicate=",Repli,sep=""), col.lab=rgb(0,0.5,0))
title(ylab=paste(names(ptdv)[pheno],":-log10(p-vlue)",sep=""),col.lab=rgb(.5,0,0))
# legend(1,5,c("H","L","tri"),cex=0.8,col=c("blue","orange","pink"),pch=c(19,17,16), lty=1:2)
dev.off()













