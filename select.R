#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
# index<-1
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
# dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
pas<- para$pas[index]
chr<-para$chr[index]
##b<-para$pheno[index]
Repli<-1000000
##pheno L position 5 in testfrm position 5 also in dcornu
b<-5
library(energy)
library(plyr)
library(foreach)
library(doMC)
##number of node to run
registerDoMC(2)
print("selectpar")
print(chr)
print(pas)
testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))
dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"dr",pas,sep=""),header=T,sep=" ")
dcornu<-cbind(c(1:nrow(dcor)),dcor)
if (!chr==4) {list<-dcornu[,b]>0.064} else {list<-dcornu[,b]>max(dcornu[,b])-.0019}
smalldcor<-dcornu[list,]
rlist<- foreach(j= 1:nrow(smalldcor), .combine='rbind') %dopar% 
{dcov.test(testfrm[,b],testfrm[, (smalldcor[j,1]+13):(smalldcor[j,1]+13+pas)],R=Repli)$p.value}
toto<-cbind(smalldcor,rlist) 
names(toto[,8])<-"pR"
outfile <- paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",pas,"select",sep="")
write.table(toto,outfile)








