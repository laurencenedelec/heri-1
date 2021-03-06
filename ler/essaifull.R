#!/usr/bin/env  Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
## dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
pas <- para$pas[index]
##pas<-15

chr<-para$chr[index]
## chr<-4
# Repli<-1500
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

snpname= read.table(paste("~/NFG/raw/Hd",chr,"fullcl.snplist",sep=""), header=F, sep='')
g<-grep("cnv",snpname[,1])
if (!length(g)==0) {snpnamesous=snpname[-g,1]} else {snpnamesous<-snpname}
 snpnamesous<-as.data.frame(snpnamesous)

testfrm=read.table(paste("~/NFG/raw/leHd",chr,"fullcl",sep=''))

 
ndonne<-(ncol(testfrm)-13)*4
ptdv=matrix(rep(0,ndonne),ncol=4)

pheno<-5

 rlist1<- foreach(j = 14:2000, .combine='rbind') %dopar% dcor.ttest(testfrm[,pheno],testfrm[,j:(j+pas)])$p.value
ptdv[1:(2000-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]  

 rlist2<- foreach(j = 20001:(ncol(testfrm)-pas), .combine='rbind') %dopar% dcor.ttest(testfrm[,pheno],testfrm[,j:(j+pas)])$p.value
ptdv[(2001-13):(ncol(testfrm)-pas-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]

snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")
ptdv<-cbind(snpnamesous,ptdv)
names(ptdv)<-c("snp","KolH","KolL","Trigl","BMI")
ptdv<-join(snpposi[,c(2,4)],ptdv, 'snp',type='right')

 write.table(ptdv,file=paste("~/NFG/result/allessaiHdap",chr,"fullpttest","pas=",pas,sep=""))













