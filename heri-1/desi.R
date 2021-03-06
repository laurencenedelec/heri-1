#!/usr/bin/env Rscript
##drawn for LDL dcor and weirdpv
index <- as.numeric(commandArgs(trailingOnly=TRUE))
# index<-1
para <- read.csv("~/NFG/ler/dessinpara.csv", header=TRUE, as.is=TRUE)
# dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
## pas <- para$pas[index]
chr<-para$chr[index]

b<-4
 
dcor=read.table("~/NFG/result/allessaiHdap4dr15",header=TRUE,sep=" ")
pdf(file= "~/NFG/result/div")
par(mfrow=c(4,1))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.1),ann=FALSE)
title(xlab="dcor 4 pas 15", col.lab=rgb(0,0.5,0))


pcvpas<-read.table("~/NFG/result/allessaiHdap19pttestpas=20",header=T,sep="")
plot(pcvpas[,2]/1000000,-log(pcvpas[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,20),ann=FALSE)
title(xlab="pas 20 weird p 19", col.lab=rgb(0,0.5,0))

pcvpas<-read.table("~/NFG/result/allessaiHdap4pttestpas=20",header=T,sep="")
plot(pcvpas[,2]/1000000,-log(pcvpas[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,20),ann=FALSE)
title(xlab="pas 20 weird p 4", col.lab=rgb(0,0.5,0))

pcvpas<-read.table("~/NFG/result/allessaiHdap1pttestpas=20",header=T,sep="")
plot(pcvpas[,2]/1000000,-log(pcvpas[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,20),ann=FALSE)
title(xlab="pas 20 weird p 1", col.lab=rgb(0,0.5,0))


dev.off()



