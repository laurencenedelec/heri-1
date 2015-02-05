#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
# index<-1
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
# dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
pas<- para$pas[index]
chr<-para$chr[index]
##b<-para$pheno[index]
Repli<-1000000
setwd("~/Documents/data/sherlock.nov")
#setwd("~/NFG")
##pheno L position 5 in testfrm position 5 also in dcornu
b<-5
library(energy)
library(plyr)
library(foreach)
library(doMC)
##number of node to run
registerDoMC(2)
#the genone data clean
testfrm=read.table(paste("raw/leHd",chr,"cl",sep=''))
#the dcor compute
dcor=read.table(paste("result/allessaiHdap",chr,"dr",pas,sep=""),header=T,sep=" ")
#the selection
dcor.number<-cbind(c(1:nrow(dcor)),dcor)
if (!chr==4) {
  list<-dcor.number[,b]>0.064
} else {
  list<-dcor.number[,b]>max(dcor.number[,b])-.0019
}
small.dcor<-dcor.number[list,]
rlist<- foreach(j= 1:nrow(small.dcor), .combine='rbind') %dopar% {
#compute boots pvalue
dcov.test(testfrm[,b],testfrm[, (small.dcor[j,1]+13):(small.dcor[j,1]+13+pas)],R=Repli)$p.value
}
res<-cbind(small.dcor,rlist) 
names(res[,8])<-"pR"
outfile <- paste("result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",pas,"select",sep="")
write.table(res,outfile)








