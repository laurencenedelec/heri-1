library(plyr)
library(energy)
library(car)
library("matrixStats")
setwd("~/NFG/ler")
source("mondcov.R")

tfam<-read.table("~/NFG/raw/NFBC_transpose.tfam")
names(tfam)<-c("SUBJID","x","y","t","z","w")

##pheno  in the same order as tfam                           
testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
#Keep only some pheno                 
pheno<-testpheno[,c(-3,-4,-5,-6,-7,-8)]                                       
pheno.order<-join(tfam[,c(1,2)],pheno[,c(2,3,4,5,6)],"SUBJID")   
rownames(pheno.order) <- paste(tfam$SUBJID)
##for those pheno add the individu id 3=K 4=H 5=L 6=T


##Keep the pheno you wants normalize pheno
pheno<-data.matrix( pheno.order[,c(4,5,6,3)])
rownames(pheno) <- paste(tfam$SUBJID)
pheno.t<-t(pheno)-colMeans(pheno)
pheno<-t(pheno.t)
pheno.sd<-colSds(pheno)
pheno.t<-t(pheno)*(1/pheno.sd)
pheno<-t(pheno.t)

#write.table(pheno,file="~/NFG/raw/pheno.normal")
#V=read.table("~/NFG/raw/pheno.normal")

source("~/NFG/ler/loadk.R")
K.divers<-load_K()

#Kdis.plink= read.table('~/NFG/result/plinkdis.mdist', sep='')
#Kdis.plink<-1-Kdis.plink



#browser()

source("~/NFG/ler/heat.R")
res<-c()

#one plot
# PlotHeatmapCorrelationMatrix( Kdis)}

##for all method, choice of K
for (j in 1:6)
{ 
#{if (j==4) {Kdis<-K.divers[[j]]} else 
Kdis<-K.divers[[j]]
id<-rownames(Kdis)
Kdis<-data.matrix(2*Kdis)

#print(colSums(is.na(Kdis)))
##detect NA in Kdis
Kdis.na<-na.omit(Kdis)
#miss<-setdiff(Kdis.na[,1],Kdis[,1])
#print("diff Kdis na ")
#print(miss)
#remove NA in Kdis
for (l in 1:ncol(Kdis)) {ok <- !is.na(Kdis[,l])
                               Kdis[!ok,l]<- mean(Kdis[ok,l])  }
Kdis<-Kdis[as.vector(id),as.vector(id)]
norm<-dcov_mc(Kdis,Kdis,1)
Kdis.list<-as.numeric(c(Kdis))


print(norm)
#print(Kdis.mean)
#for all pheno 
for (i in 1:4) 
{print(i)
  Y<-pheno[,i]
  names(Y) <- paste(tfam$SUBJID)
  
## filter the pheno with the K$id
#if (!identical(id,tfam$SUBJID)) { 
Y<-Y[as.vector(id)]
#} 
Y.produit<-Y %*%t(Y)
Y.produit.list<-as.numeric(c(data.matrix(Y.produit)))
nle<-length(Y)
dia<-c()
for (k in 1:nle){ dia<-cbind(dia,c((k-1)*nle+k)) }
dia<-c(dia)
Y.produit.sans<-Y.produit.list[-dia]
##compute the heritability-linear regression
Kdis.list.sans<-Kdis.list[-dia]
Kdis.list.sans<-as.numeric(as.vector(Kdis.list.sans ))
Y.produit.list.sans<-as.numeric(as.vector(  Y.produit.sans))
heri.lm.sans<-lm(Y.produit.list.sans~Kdis.list.sans)
print(length(Kdis.list.sans))
heri.sans<-summary(heri.lm.sans)$coef[2,1]
pv.heri.sans<-summary(heri.lm.sans)$coef[2,4]
Kdis.mean<-mean(Kdis.list.sans)
if (j==1)    red<- Kdis.list.sans<(1.2)*Kdis.mean  
if (j==2)    red<- Kdis.list.sans<(1.2)*(1-Kdis.mean)
if (j==4)    red<- Kdis.list.sans<(1.2)*(1-Kdis.mean)
if (j==3)    red<- Kdis.list.sans<(1.2)*(1-Kdis.mean)
if (j==5)    red<- Kdis.list.sans<(1.5)*Kdis.mean
if (j==6)    red<- Kdis.list.sans<(1.5)*Kdis.mean

Kdis.list.red<-Kdis.list.sans[red]
Y.produit.list.red<- Y.produit.list.sans[red]
heri.lm.red<-lm(Y.produit.list.red~Kdis.list.red)


heri.lm.all<-lm(Y.produit.list~Kdis.list)
heri.red<-summary(heri.lm.red)$coef[2,1]
pv.heri.red<-summary(heri.lm.red)$coef[2,4]
#heri.red<-1
#pv.heri.red<-1

heri.all<-summary(heri.lm.all)$coef[2,1]
pv.heri.all<-summary(heri.lm.all)$coef[2,4]

##compute the dcov 
#ny<-length(Y)
#Ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Yp

heri.dcov<-dcov_mc(Y.produit,Kdis,1)
heri.dcov<-heri.dcov/norm

#ne marche pas
#norm.dcov.energy<-(dcov(as.dist(Kdis),as.dist(Kdis),1))^2
#heri.dcov.energy<-(dcov(as.dist(Y.produit),as.dist(Kdis),1))^2
#heri.dcov.energy<-heri.dcov.energy/norm.dcov.energy

#heri.dcov.energy<-1
#norm.dcov<-1

##the results
herita<-c(j,i,heri.sans,pv.heri.sans,heri.red,pv.heri.red,heri.all,pv.heri.all,heri.dcov, Kdis.mean,)
res<-rbind(res,herita)


}
}
res<-data.matrix(res)
colnames(res)<-c('method', 'pheno','sansdiag','pvsansdiag', 'lm.red','pv.lm.red', 'lm','pv.lm','heri.dcov','meanK')
rownames(res)<-c('K.ibs','K.idb','K.ibs','K.cgta',"k ","kk",'H.ibs','H.idb','H.ibs','H.cgta',"h",'hh",'L.ibs','L.idb','L.ibs','L.cgta',"L",'ll','T.ibs','T.idb','T.ibs','T.cgta','T','tt')
write.table(res,file="~/NFG/ler/res5.heri.RData", row.names=TRUE, col.names=TRUE)



