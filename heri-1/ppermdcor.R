#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
pas <- para$pas[index]
chr<-para$chr[index]
Repli<-1000

testped = read.table(paste("~/NFG/raw/Hd",chr,"clA.raw",sep=""), header=T,sep='')
testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
testfi<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))

library(energy)
library(plyr)

## testfn<-testfi
## for (j in 2:ncol(testfi)){is.na(testfn[,j]) <- testfn[,j]=="X"}
## testfnr<-testfn[,c(-15,-17,-18,-31)]
## for (j in 30:ncol(testfnr)) {ok <- !is.na(testfnr[,j]) 
##     			       testfnr[!ok,j]<- mean(testfnr[ok,j])  } 
## colSums(is.na(testfnr)) 
## testfrm<-testfnr[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
## testfrm<-na.omit(testfrm)

snpname= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
 g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)
## cut<-( (nrow(snpnamesous)-9):nrow(snpnamesous))
## snpsoussous<-snpnamesous[-cut,]
## snpsoussous<-as.data.frame(snpsoussous)
## write.table(testfrm,file=paste("~/NFG/raw/leHd",chr,"cl",sep=""))

testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))

## gfrm<-grep("cnv",names(testfrm))
## testfrm<-testfrm[,-gfrm]
 


dnul<-rep(1==0,ncol(testfrm))
for (j in 13:ncol(testfrm)) 
{ dnul[j]<-max(testfrm[ ,j])==0}
bnul<-dnul[-(1:13)]
testfrm<-testfrm[,!dnul]
snpnamesous<-snpnamesous[!bnul,]
snpnamesous<-as.data.frame(snpnamesous)

ndonne<-(ncol(testfrm)-13)*4
dc=matrix(rep(0,ndonne),ncol=4)
ptdv=matrix(rep(0,ndonne),ncol=4)

for (j in 14:(ncol(testfrm)-pas) )
     {  	
     tempcor<-dcov.test(testfrm[,4],testfrm[,j:(j+pas)],R=Repli)
     ptdv[j-13,1]<-tempcor$p.value
     }

for (j in 14:(ncol(testfrm)-pas) )
     {  	
     tempcor<-dcov.test(testfrm[,5],testfrm[,j:(j+pas)],R=Repli)
     ptdv[j-13,2]<-tempcor$p.value
     }

for (j in 14:(ncol(testfrm)-pas) )
     {  	
     tempcor<-dcov.test(testfrm[,6],testfrm[,j:(j+pas)],R=Repli)
     ptdv[j-13,3]<-tempcor$p.value
     }
 
snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")


ptdv<-cbind(snpnamesous,ptdv)
dc<-cbind(snpnamesous,dc)

names(ptdv)<-c("snp","KolH","KolL","Trigl","BMI")
names(dc)<-c("snp","KolH","KolL","Trigl","BMI")
ptdv<-join(snpposi[,c(2,4)],ptdv, 'snp',type='right')
dc<-join(snpposi[,c(2,4)],dc,'snp', type='right')


write.table(ptdv,file=paste("~/NFG/result/Hdap",chr,"pperm",pas,sep=""))
write.table(dc,file=paste("~/NFG/result/Hdap",chr,"dcper",pas,sep=""))

## dessiner les pvalues
for (pheno in 3:6 )
{
jpeg(file=paste("~/NFG/result/Hd",pheno,"-",chr,"pperm",pas,".jpg",sep=""))
plot(ptdv[,2]/1000000,-log(ptdv[,pheno],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,13),ann=FALSE)
title(xlab=paste("Chr",chr,"position(Mb) by dc-per pas=",pas,"pheno=",pheno,"replicate=",Repli), col.lab=rgb(0,0.5,0) )
title(ylab=paste(names(ptdv)[pheno],":-log10(p-vlue)",sep=""),col.lab=rgb(.5,0,0))
dev.off()
}













