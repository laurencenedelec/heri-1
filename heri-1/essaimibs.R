#!/usr/bin/env  Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
## dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")
##pas <- para$pas[index]
##pas<-15
##pas<-20
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
##registerDoMC(20)

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
# rlist<- foreach(j = 14:(ncol(testfrm)-pas), .combine='rbind') %dopar% dcor.ttest(testfrm[,pheno],testfrm[,j:(j+pas)])$p.value
#ptdv[1:(ncol(testfrm)-pas-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]  

## list<-foreach(j=14:(ncol(testfrm)-pas),.combine='rbind') %dopar% pcov(testfrm[,pheno],testfrm[,j]) 
## pcv[1:(ncol(testfrm)-pas-13),ipheno] <- list[1:(ncol(testfrm)-pas-13),1]

##list<-foreach(j=14:(ncol(testfrm)-pas),.combine='rbind') %dopar% dcor(testfrm[,pheno],testfrm[,j:(j+pas)],1)
##dc[1:(ncol(testfrm)-pas-13),ipheno] <- list[1:(ncol(testfrm)-pas-13),1]


#rlist<- foreach(j = 14:(ncol(testfrm)-pas),.combine='rbind') %dopar% pcov_vvdim(testfrm[,pheno],testfrm[,j:(j+pas)]) 
#pcvpas[1:(ncol(testfrm)-pas-13),ipheno] <- rlist[1:(ncol(testfrm)-pas-13),1]

 mibs=read.table('~/NFG/ler/plink.mibs', header=F , sep='')
 dmibs<-1-mibs
if (chr==19) {missindi<-c(1737,4756)} 
if (chr==4) {missindi<-c(158,171,452,869,1030,1466,1584,1682,2007, 2230, 2476, 3445, 4121, 4354, 4617,4642, 4858, 4930, 5179, 5214, 5332, 5346)}
if (chr==1)
{missindi<-c(62,64,102,109,209,227,288,398,435,507,529,536,567,575,590,594,618,672,769,826,837,900,975,986, 990,1040,1064,1128,1133,1154,1192,1260,1274,1287,1321)
missindi<-c(missindi,1332,1333,1341,1342,1348,1395,1405,1448,1462,1499,1560,1680,1751,1775,1779,1790,1860,1889,1890,1965,2045,2176,2195,2200, 2282,2297, 2388, 2407)
missindi<-c(missindi, 2466, 2482, 2587, 2590,2595, 2653,2657, 2694, 2710, 2715, 2779, 2793,2865, 2884, 2893, 2931, 2945, 2952, 3025, 3034, 3059, 3208, 3217, 3258)
missindi<-c(missindi, 3274, 3289, 3328,3384, 3401, 3537,3538,3587,3626,3650,3697,3728,3747,3749,3783,3819,3906,3986,4007, 4016, 4063, 4084, 4114, 4160, 4191, 4215)
missindi<-c(missindi, 4407, 4430, 4432, 4489, 4518,4549,4607,4635, 4668, 4687, 4696, 4714, 4764, 4826, 4924, 4993, 5100, 5118, 5121, 5133, 5186, 5248, 5294, 5307)
missindi<-c(missindi, 5327, 5339, 5355, 5374, 5396)
}

out<- dcor_vm(testfrm[,pheno],dmibs[-missindi, -missindi],1)

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
# write.table(pcv,file=paste("~/NFG/result/allessaiHdap",chr,"pcvsecond",pas,sep=""))
# write.table(pcvpas,file=paste("~/NFG/result/allessaiHdap",chr,"pdrmoi",pas,sep=""))
##write.table(dc,file=paste("~/NFG/result/allessaiHdap",chr,"dr",pas,sep=""))
write.table(out,file=paste("~/NFG/result/allessaiHdap",chr,"drmibs",sep=""))
#write.table(dcallele,file=paste("~/NFG/result/allessaiHdap",chr,"drallele",sep=""))
# write.table(ptdv,file=paste("~/NFG/result/allessaiHdap",chr,"pttest","pas=",pas,sep=""))













