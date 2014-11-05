#!/usr/bin/env Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
# index<-1
para <- read.csv("~/NFG/ler/dessinpara.csv", header=TRUE, as.is=TRUE)
# dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
## pas <- para$pas[index]
chr<-para$chr[index]
##b<-para$pheno[index]
Repli<-1000000
pas<-0
b<-4
##pvalue=read.table(paste("~/NFG/result/Hdap",chr,"pv",pas,sep=""),header=T,sep=" ")
##jpeg(file= paste("~/NFG/result/H",b,"-ch",chr,"pv",pas,".jpg",sep=""))
##plot(pvalue[,2]/1000000,-log(pvalue[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,18),ann=FALSE)
##title(xlab=paste("Chr",chr,"pas=",pas,"pheno=",b,"position (MB)by dcov",sep=""), col.lab=rgb(0,0.5,0))
##title(ylab="-log10(p-vlue)",col.lab=rgb(.5,0,0)) 
### legend(1,5,c("H","L","tri"),cex=0.8,col=c("blue","orange","pink"),pch=c(19,17,16), lty=1:2)
##dev.off()
##dcor=read.table(paste("~/NFG/result/Hdap",chr,"dc",pas,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
##plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.2),ann=FALSE)
##title(xlab=paste("Chr",chr,"pas=",pas,"pheno=",names(dcor)[b],"position (MB)by dcov",sep=""), col.lab=rgb(0,0.5,0))
##title(ylab="d-cor",col.lab=rgb(.5,0,0))
##dev.off()
##ptdv=read.table(paste("~/NFG/result/allessaiHdap",chr,"pR",pas,sep=""))
##pdf(file=paste("~/NFG/result/allessaiHd",names(ptdv)[b],"-",chr,"pR",pas,".pdf",sep=""))
##plot(ptdv[,2]/1000000,-log(ptdv[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,5),ann=FALSE)
##title(xlab=paste("Chr",chr, "position(Mb) by dcov-pas=",pas,"pheno=",names(ptdv)[b],"repli=",500,sep=""), col.lab=rgb(0,0.5,0))
##title(ylab=paste(names(ptdv)[b],":-log10(p-vlue)",sep=""),col.lab=rgb(.5,0,0))
##dev.off()

pvalue=read.table(paste("~/NFG/resultv/Hdap",chr,"pstand",sep=""),header=T,sep=" ")
pdf(file= paste("~/NFG/result/toutH1", names(pvalue)[b],"-ch",chr,".pdf",sep=""))
#par(mfrow=c(3,2), mar=c(1.1,4.1,3.5,2.1))
#if (chr==4) {par(mfrow=c(4,1))} else {par(mfrow=c(3,1))}
par(mfrow=c(3,1))
plot(pvalue[,2]/1000000,-log(pvalue[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,10),ann=FALSE)
sel5=read.table(paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",5,"select",sep=""))
points(sel5[,3]/1000000,-log(sel5[,8],10) ,col="blue",pch=19)
if (!chr==4) {sel10=read.table(paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",10,"select",sep=""))
points(sel10[,3]/1000000,-log(sel10[,8],10) ,col="green",pch=18)   }
title(xlab=paste("Chr",chr,"pheno=",b,"position (MB)by stand and p-repli-dcor 5 10",sep=""), col.lab=rgb(0,0.5,0))
title(ylab=paste(names(pvalue)[b],":-log10(p-vlue) stand",sep=""),col.lab=rgb(.5,0,0))

##sel5=read.table(paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",5,sep=""))
##sel10=read.table(paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",5,sep=""))
#### pdf(file=paste("~/NFG/result/allessaiHd",names(ptdv)[b],"-",chr,"pR","repli=",Repli,"pas=",pas,".pdf",sep=""))
##plot(sel5[,2]/1000000,-log(sel[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,5),ann=FALSE)
##plot(sel10[,2]/1000000,-log(sel[,b],10) ,col="green",pch=18,xlim=c(40,90),ylim=c(0,5),ann=FALSE)
##title(xlab=paste("Chr=",chr, " position(Mb) by dcov-pas=",pas," pheno=",names(sel)[b]," repli=",500,sep=""), col.l##ab=rgb(0,0.5,0))
##title(ylab=paste(names(sel)[b],":-log10(p-vlue) repli",sep=""),col.lab=rgb(.5,0,0))

##pcvpas<-read.table(paste("~/NFG/result/allessaiHdap",chr,"pcvsecond",pas,sep=""))
##pdf(file=paste("~/NFG/result/allessaiHd",names(pcvpas)[b],"-",chr,"pcvsecond",pas,".pdf",sep=""))
##plot(pcvpas[,2]/1000000,-log(pcvpas[,b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,4),ann=FALSE)
##title(xlab=paste("Chr",chr, "position(Mb) by dcov-pas=",pas,"pheno=",names(pcvpas)[b]," seconddcov",sep=""), col.l##ab=rgb(0,0.5,0))
##title(ylab=paste(names(pcvpas)[b],"(p-vlue) second",sep=""),col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"dr",5,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.1),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=",5," pheno=",names(dcor)[b]," position (MB) dcor",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dcor pas 5",col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"dr",10,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.1),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=",10," pheno=",names(dcor)[b]," position (MB) dcor",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dcor pas 10",col.lab=rgb(.5,0,0))



dev.off()

pdf(file= paste("~/NFG/result/toutH2", names(pvalue)[b],"-ch",chr,".pdf",sep=""))
#par(mfrow=c(3,2), mar=c(1.1,4.1,3.5,2.1))
if (chr==4) {par(mfrow=c(3,1))} else {par(mfrow=c(2,1))}
 if (chr==4) {
dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"dr",15,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.1),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=",10," pheno=",names(dcor)[b]," position (MB) dcor",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dcor pas 15",col.lab=rgb(.5,0,0))
}
dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"pttestpas=",20,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) pttest pas 20",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="pttest",col.lab=rgb(.5,0,0))


dcor=read.table(paste("~/NFG/result/Hdap",chr,"pv",10,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) pttest pas 10",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="pttest",col.lab=rgb(.5,0,0))

dev.off()

if (chr==19) {
pdf(file= paste("~/NFG/result/toutH4", names(pvalue)[b],"-ch",chr,".pdf",sep=""))
#par(mfrow=c(3,2), mar=c(1.1,4.1,3.5,2.1))
if (chr==4) {par(mfrow=c(3,1))} else {par(mfrow=c(2,1))}
 if (chr==4) {
dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"dr",15,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(0,90),ylim=c(0,.1),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=",10," pheno=",names(dcor)[b]," position (MB) dcor",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dcor pas 15",col.lab=rgb(.5,0,0))
}
dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"fullpttestpas=",20,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(0,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) fullpttest pas 20",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="pttest",col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"fullpttestpas=",10,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(0,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) fullpttest pas 10",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="pttest",col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"fullpttestpas=",5,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(0,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) fullpttest pas 5",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="pttest",col.lab=rgb(.5,0,0))

dev.off()
}

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"pdrmoi",10,sep=""),header=T,sep=" ")

pdf(file= paste("~/NFG/result/toutH3",names(dcor)[b],"-ch",chr,"pas=10",".pdf",sep=""))
par(mfrow=c(3,1))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=10"," pheno=",names(dcor)[b]," position (MB) pmoi pas 10",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="psol",col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"pdrmoi",20,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,-log(dcor[, b],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,10),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) pmoi pas 20",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="psol",col.lab=rgb(.5,0,0))

dcor=read.table(paste("~/NFG/result/allessaiHdap",chr,"drmoi",10,sep=""),header=T,sep=" ")
##pdf(file= paste("~/NFG/result/allH",names(dcor)[b],"-ch",chr,"dc",pas,".pdf",sep=""))
plot(dcor[,2]/1000000,dcor[, b] ,col="red",pch=18,xlim=c(40,90),ylim=c(0,.1),ann=FALSE)
title(xlab=paste("Chr=",chr," pas=20"," pheno=",names(dcor)[b]," position (MB) drmoi pas 10",sep=""), col.lab=rgb(0,0.5,0))
title(ylab="dr",col.lab=rgb(.5,0,0))

dev.off()


