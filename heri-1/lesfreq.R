library(plyr)
library(energy)
library(car)
library("matrixStats")
print("lesfreq")
gene<-"Hd1clA"
gelis<-"Hd1cl"
##gene<-"allgene"
testped = read.table(paste('~/NFG/raw/',gene,'.raw',sep=""), header=T,sep='')
testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
 testfn<-merge(testpheno,testped, by.x=c('SUBJID'),by.y=c('IID'))

####remove NA by mean on the genome
for (j in 2:ncol(testfn)){is.na(testfn[,j]) <- testfn[,j]=="X"}
testfn<-testfn[,c(-15,-17,-18,-31)]
for (j in 30:ncol(testfn)) {ok <- !is.na(testfn[,j])
    testfn[!ok,j]<- mean(testfn[ok,j])  }

testfn<-testfn[,c(-4,-5,-6,-7,-8,-9,-13,-14,-15,-16,-17,-18,-19,-23,-24)]
test<-na.omit(testfn)
##normalize the geno
tutu<-test[,14:ncol(test)]
tutu<-data.matrix(tutu)
tutu[,ncol(tutu)]<-tutu[,ncol(tutu)]-colMeans(tutu)[ncol(tutu)]
test<-cbind(test[1:13],tutu)  


write.table(test,file=paste("~/NFG/raw/leall",gene,"gene",sep=""))
## test=read.table('~/NFG/raw/leallgene',sep='')
## remove cnv for data
gtest<-grep("cnv",names(test))
if (!length(gtest)==0) {test=test[,-gtest]} else {test<-test}


## Mettre les nom au bon format enlever des -A autre facon ?
snpnames= read.table( paste("~/NFG/raw/",gelis,".snplist",sep=""), header=F, sep='')
##remove cnv from name
g<-grep("cnv",snpnames[,1])
if (!length(g)==0) {snpnames=snpnames[-g,1]} else {snpnames<-snpnames}
 snpnames<-as.data.frame(snpnames)


##snpnames<-names(test)
##snpnames<-sub('_G$', '', snpnames[])
##snpnames<-sub('_A$', '', snpnames[])
##snpnames<-sub('_0$', '', snpnames[])
##snpnames<-sub('_C$', '', snpnames[])
##snpnames=snpnames[-(1:13)]
#### snpnames<-as.data.frame(snpnames)

## remove zero from data and name
 dnul<-rep(1==0,ncol(test))
 for (j in 13:ncol(test))
 { dnul[j]<-max(test[ ,j])==0}
 dnul<-as.logical(dnul)
 bnul<-dnul[-(1:13)]
 test<-test[,!dnul]
 snpnames<-snpnames[!bnul,]
 snpnames<-as.data.frame(snpnames)

##compute the frequency, write them
frequ=read.table('~/NFG/result/freqNFG.frq',header=T, sep='')
names(frequ)<-c("chr","snp","A1","A2","MAF","NCHROBS")

##or another way to get the frequency
freqdo<-colMeans(test[,-(1:13)])/ 2
 un<-rep(1,ncol(test)-13)
frequbis<-un-freqdo
frequbis<-as.matrix(frequbis)
##frequbis<-t(frequbis)
##yet another way qui ne marche pas
fretrois<-un
for (j in 1:(ncol(test)-13))
 {zero<-(test[,j+13]==0)
 une<-(test[,j+13]==1)
 deux<-(test[,j+13]==2)
rzero<-length(test[zero,j+13])
rune<- length(test[une,j+13])
rdeux<- length(test[deux,j+13])
 fretrois[j]<-1-(rzero+rune/2 )/(rzero+rune+rdeux )
 }

fretrois<-as.matrix(fretrois)

write.table(fretrois,file=paste("~/NFG/result/fretrois",gene,sep=""))
## add the position to the two frequency in the order of test
fre<-cbind(snpnames,frequbis,fretrois)
names(fre)<-c("snp","frequbis","fretrois")
fre<-join(fre,frequ[,c(2,5)],"snp",type='right')
snpposi=read.table('~/NFG/raw/snpposi.bim', header=F , sep='')
names(snpposi)<-c("x","snp","y","posi","z","w")
fre<-join(fre,snpposi[,c(2,4)], 'snp',type='right')
names(fre)<-c("snp","freqdo","fretrois","frequ","posi")
write.table(fre,file=paste("~/NFG/result/freq",gene,sep=""))
##analyse-draw

pdf(file="~/NFG/result/frequence.pdf")
 plot(fre[,5]/10000,fre[,2] ,col="red",pch=18,xlim=c(0,10000),ylim=c(0,3),ann=FALSE)
 points( fre[,5]/10000,1+fre[,3],col="blue",pch=19)
points( fre[,5]/10000,2+fre[,4],col="green",pch=20)
 title(xlab="snp", col.lab=rgb(0,0.5,0))
title(ylab="frequence",col.lab=rgb(.5,0,0))
 dev.off()