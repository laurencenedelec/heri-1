#!/usr/bin/env  Rscript
index <- as.numeric(commandArgs(trailingOnly=TRUE))
para <- read.csv("~/NFG/ler/para.csv", header=TRUE, as.is=TRUE)
## dataFile <- paste("~/NFG/ler/raw", para$dataFile[index], sep="/")

# parameter

# pas <- para$pas[index]
# pas<-15
# pas<-20
pas<-10
chr<-para$chr[index]
# chr<-4
 Repli<-1500
# 4,5,6,9 are K L T BMI
num.pheno<-5
# path
setwd("~/NFG/ler")

# Load libraries and source helper functions --------------------------------------
source("mondcov.R")

library(energy)
library(plyr)
library(foreach)

# Setup for parallel loop with 20 cores
# number of node to run
library(doMC)
registerDoMC(20)


# dowload the four files : pheno geno and snp.name snp.posi
snp.name= read.table(paste("~/NFG/raw/Hd",chr,"cl.snplist",sep=""), header=F, sep='')
snp.posi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snp.posi)<-c("x","snp","y","posi","z","w")

# clean cnv form snp
g<-grep("cnv",snp.name[,1])
if (!length(g)==0) {snp.name.sous=snp.name[-g,1]
} else {snp.name.sous<-snp.name}
snp.name<-as.data.frame(snp.name.sous)


geno= read.table(paste("~/NFG/raw/Hd",chr,"clA.raw",sep=""), header=T,sep='')
pheno= read.table('~/NFG/raw/Pheno', header=TRUE, sep='',skip=11)

gsub("X","NA",pheno)
pheno<-data.matrix(pheno)
# selection of less pheno
less<-c(-4,-5,-6,-7,-8,-9,-13,-14,-15, -16,-17,-18,-19,-20,-21,-22,-26,-27,-31)
pheno<-pheno[,less]

# 4,5,6,9 are K L T BMI
data<-merge(pheno,geno, by.x=c('SUBJID'),by.y=c('IID'))
# the geno start at the 14 columm
# cleanning the NA in pheno  and geno ( instead of "X")
data.clean<-data
for (j in 2:ncol(data)){is.na(data.clean[,j]) <- data.clean[,j]=="X"}
# cleanning of the geno part NA become mean
for (j in 30:ncol(data.clean)) {ok <- !is.na(data.clean[,j])
    data.clean[!ok,j]<- mean(data.clean[ok,j])  }

# for testing 0 is the correct answer
# colSums(is.na(data.clean))


## Selection that involve lost of individu rm the NA in pheno

# To keep track of the lost one
data.clean<-na.omit(data)
miss<-setdiff(data[,1],data.clean[,1])
write.table(miss,paste('~/NFG/raw/missindividu',chr,sep="") )
data<-data.clean



# testfrm=read.table(paste("~/NFG/raw/leHd",chr,"cl",sep=''))


# clean cnv from data
grep.data<-grep("cnv",names(data))
if (!length(grep.data)==0) {data<-data[,-grep.data]}


# compare size in case
ncol(data)
nrow(snp.name)

# cleaning the zero column
# dnul<-rep(1==0,ncol(data))
# for (j in 13:ncol(data))
# { dnul[j]<-max(data[ ,j])==0}
# bnul<-dnul[-(1:13)]
# if (!length(dnul)==0) {data=data[,!dnul]}
# if (!length(bnul)==0) {snp.name= snp.name[!bnul,]}
# snp.name<-as.data.frame(snp.name)






 rlist<- foreach(j=14:(ncol(data)-pas), .combine='rbind') %dopar%
# choose the object to be compute

# dcor.ttest(data[,num.pheno],data[,j:(j+pas)])$p.value
#pcov(data[,num.pheno],data[,j])
#dcor(data[,num.pheno],data[,j:(j+pas)],1)
#dcor_vvdim(data[,pheno],(data[,j:(j+pas)]))

res <- rlist[,1]




# to do
# mibs=read.table('~/NFG/ler/plink.mibs', header=F , sep='')
# dmibs<-1-mibs
# out<- dcor-vm(data[,pheno],dmibs[-c(1737,4756), -c(1737,4756)],1)
# dcallele[1:(ncol(data)-pas-13),index.pheno] <- list[1:(ncol(data)-pas-13),1]



# The pas create a blank at the end
cut<-( (nrow(snp.name)-pas+1):nrow(snp.name))
snp.name.pas<-snp.name[-cut,]

res<-cbind(snp.name.pas,res)
names(res)<-c("snp",names(pheno)[num.pheno])


res<-join(snp.posi[,c(2,4)],res, 'snp',type='right')

# choose where to write it

# write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",pas,sep=""))
# write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"pcvsecond",pas,sep=""))
  write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"-",names(pheno)[num.pheno],"drmoi",pas,sep=""))
# write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"dr",pas,sep=""))
# write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"drallele",sep=""))
# write.table(res,file=paste("~/NFG/result/allessaiHdap",chr,"pttest","pas=",pas,sep=""))













