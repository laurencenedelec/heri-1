#!/usr/bin/env Rscript
 index <- as.numeric(commandArgs(trailingOnly=TRUE))
# index<-1
para <- read.csv("~/NFG/ler/dessinpara.csv", header=TRUE, as.is=TRUE)
# dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
pas <- para$pas[index]
chr<-para$chr[index]
b<-para$pheno[index]
for (b in 3:6)
{pvalue=read.table(paste("~/NFG/result/Hdap",chr,"dc",pas,sep=""),header=T,sep=" ")
pdf(file= paste("~/NFG/result/H",b,"-ch",chr,"pv",pas,".pdf",sep=""))
plot(pvalue[,2]/1000000,pvalue[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,18),ann=FALSE)
# points(-log(pvalue[,3],10),col="blue",pch=19)
# points(-log(pvalue[,5],10),col="orange",pch=17)
# points(-log(pvalue[,6],10),col="pink",pch=16)
title(xlab=paste("Chr",chr,"pas=",pas,"pheno=",b,"position (MB)by dcov",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dcor",col.lab=rgb(.5,0,0))
# legend(1,5,c("H","L","tri"),cex=0.8,col=c("blue","orange","pink"),pch=c(19,17,16), lty=1:2)
dev.off()


pvaluesin=read.table(paste("~/NFG/result/Hdap",chr,"dc-sin",pas,sep=""),header=T,sep=" ")
jpeg(file= paste("~/NFG/result/H",b,"-ch",chr,"dcp-sin",".jpg",sep=""))
plot(pvaluesin[,2]/1000000,-log(pvaluesin[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,18),ann=FALSE)
# points(-log(pvalue[,3],10),col="blue",pch=19)
# points(-log(pvalue[,5],10),col="orange",pch=17)
# points(-log(pvalue[,6],10),col="pink",pch=16)
title(xlab=paste("Chr",chr,"pas=",pas,"pheno=",b,"position (MB)by dcov",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="-log10(p-vlue)",col.lab=rgb(.5,0,0))
# legend(1,5,c("H","L","tri"),cex=0.8,col=c("blue","orange","pink"),pch=c(19,17,16), lty=1:2)
dev.off()

}