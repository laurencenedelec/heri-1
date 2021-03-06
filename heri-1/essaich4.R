#!/usr/bin/env  Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
## dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
##pas <- para$pas[index]
pas<-15
##chr<-para$chr[index]
 chr<-4
 Repli<-1500
setwd("~/NFG/ler")

# Load libraries and source helper functions --------------------------------------
source("mondcov.R")


library(energy)
library(plyr)
print("essai")

# Setup for parallel loop with 20 cores
library(foreach)
library(doMC)
##number of node to run
registerDoMC(20)

snpname= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)

testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))

##cleaning the zero column
##dnul<-rep(1==0,ncol(testfrm))
##for (j in 13:ncol(testfrm)) 
##{ dnul[j]<-max(testfrm[ ,j])==0}
##bnul<-dnul[-(1:13)]
##if (!length(dnul)==0) {testfrm=testfrm[,!dnul]} else {testfrm<-testfrm}
##if (!length(bnul)==0) {snpnamesous= snpnamesous[!bnul,]} else {snpnamesous<-snpnamesous}
## snpnamesous<-as.data.frame(snpnamesous)

 
ndonne<-(ncol(testfrm)-13)*4
ptdv=matrix(rep(0,ndonne),ncol=4)
pcv=matrix(rep(0,ndonne),ncol=4)
pcvpas=matrix(rep(0,ndonne),ncol=4)
dc=matrix(rep(0,ndonne),ncol=4)
dcallele=matrix(rep(0,ndonne),ncol=4)
ipheno<-0  
##for (pheno in c(4,5,6,9))
##{ipheno<-ipheno+1
pheno<-5
ipheno<-2
## rlist<- foreach(j = 14:(ncol(testfrm)-pas), .combine='rbind') %dopar% dcov.test(testfrm[,pheno],testfrm[,j:(j+pas)],R=Repli)$p.value
## ptdv[1:(ncol(testfrm)-pas-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]  

## list<-foreach(j=14:(ncol(testfrm)-pas),.combine='rbind') %dopar% pcov(testfrm[,pheno],testfrm[,j]) 
## pcv[1:(ncol(testfrm)-pas-13),ipheno] <- list[1:(ncol(testfrm)-pas-13),1]

list<-foreach(j=14:(ncol(testfrm)-pas),.combine='rbind') %dopar% dcor(testfrm[,pheno],testfrm[,j:(j+pas)],1)
dc[1:(ncol(testfrm)-pas-13),ipheno] <- list[1:(ncol(testfrm)-pas-13),1]


##rlist<- foreach(j = 14:(ncol(testfrm)-pas),.combine='rbind') %dopar% corperdim(testfrm[,pheno],testfrm[,j:(j+pas)]) 
##pcvpas[1:(ncol(testfrm)-pas-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]

##mibs=read.table('~/NFG/ler/plink.mibs', header=F , sep='')
##dmibs<-1-mibs
##list<-foreach(j=14:(ncol(testfrm)-pas),.combine='rbind') %dopar% dcorvm(testfrm[,pheno],dmibs,1)
##dcallele[1:(ncol(testfrm)-pas-13),ipheno] <- list[1:(ncol(testfrm)-pas-13),1]  

##}

snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")

ptdv<-cbind(snpnamesous,ptdv)
pcv<-cbind(snpnamesous,pcv)
dc<-cbind(snpnamesous,dc)
pcvpas<-cbind(snpnamesous,pcvpas)
dcallele<-cbind(snpnamesous,dcallele)
names(dcallele)<-c("snp","KolH","KolL","Trigl","BMI")
names(dc)<-c("snp","KolH","KolL","Trigl","BMI")
names(pcvpas)<-c("snp","KolH","KolL","Trigl","BMI")
names(ptdv)<-c("snp","KolH","KolL","Trigl","BMI")
names(pcv)<-c("snp","KolH","KolL","Trigl","BMI")

dc<-join(snpposi[,c(2,4)],dc, 'snp',type='right')
ptdv<-join(snpposi[,c(2,4)],ptdv, 'snp',type='right')
pcv<-join(snpposi[,c(2,4)],pcv,'snp', type='right')
pcvpas<-join(snpposi[,c(2,4)],pcvpas,'snp', type='right')
dcallele<-join(snpposi[,c(2,4)],dc, 'snp',type='right')
## write.table(ptdv,file=paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",pas,sep=""))
## write.table(pcv,file=paste("~/NFG/result/allessaiHdap",chr,"pcvsecond",pas,sep=""))
## write.table(pcvpas,file=paste("~/NFG/result/allessaiHdap",chr,"drmoi",pas,sep=""))
write.table(dc,file=paste("~/NFG/result/allessaiHdap",chr,"dr",pas,sep=""))
##write.table(dcallele,file=paste("~/NFG/result/allessaiHdap",chr,"drallele",sep=""))














