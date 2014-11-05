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
##for those pheno add the individu id 3=K 4=H 5=L 6=T
Kpheno<-data.matrix(pheno.order[,c(1,2,3)])

#write.table(Kpheno,file="~/NFG/raw/Kpheno",row.names=F,col.names=F)
Hpheno<-data.matrix(pheno.order[,c(1,2,4)])                                        
# write.table(Hpheno,file="~/NFG/raw/Hpheno",row.names=F,col.names=F)
Lpheno<-data.matrix(pheno.order[,c(1,2,5)])                                        
# write.table(Lpheno,file="~/NFG/raw/Lpheno",row.names=F,col.names=F)
Tpheno<-data.matrix(pheno.order[,c(1,2,6)])
# write.table(Tpheno,file="~/NFG/raw/Tpheno",row.names=F,col.names=F)    

# Hpheno<- read.table("~/NFG/raw/Hpheno")
# Lpheno<-read.table("~/NFG/raw/Lpheno")
# Tpheno<-read.table("~/NFG/raw/Tpheno")
# Hpheno<-data.matrix(Hpheno)
# Lpheno<-data.matrix(Lpheno)
# Tpheno<-data.matrix(Tpheno)

##Keep the pheno you wants normalize pheno
pheno<-matrix(c(Hpheno[,3],Lpheno[,3],Tpheno[,3],Kpheno[,3]),ncol=4)
pheno.t<-t(pheno)-colMeans(pheno)
pheno<-t(pheno.t)
pheno.sd<-colSds(pheno)
pheno.t<-t(pheno)*(1/pheno.sd)
pheno<-t(pheno.t)
#write.table(pheno,file="~/NFG/raw/pheno.normal")
#V=read.table("~/NFG/raw/pheno.normal")

#source("~/NFG/ler/loadk.R")
#K.divers<-load_K()
K.ibsg<-read.table("~/NFG/raw/K.ibsg")
K.ibsg.order<-
K.ibdg<-read.table("~/NFG/raw/K.ibdg")
K.ibdg.order<-
Kdis.plink= read.table('~/NFG/result/plinkdis.mdist', sep='')
Kdis.plink<-1-Kdis.plink
Kdis.plink.order<-
K.gcta<-read.table("~/NFG/raw/kpgc")
K.gcta.order<-
#all method together
K.divers<-list(K.ibsg,K.ibdg,Kdis.plink,K.gcta)
print("order id in CGTA")
identical(K.gcta$id$V1,tfam$SUBJID)

#browser()
source("~/NFG/ler/heat.R")

#one plot
#for (j in 1:4)
#{if (j==4) {Kdis<-data.matrix(K.divers[[j]])} else {Kdis<-data.matrix( 2*K.divers[[j]] )}
# PlotHeatmapCorrelationMatrix( Kdis)}


res<-c()

##for all method
for (j in 1:4)
{ 
#{if (j==4) {Kdis<-data.matrix(K.divers[[j]])} else 
{Kdis<-data.matrix( 2*K.divers[[j]])} 
print('j=')
print(j)
#print(colSums(is.na(Kdis)))
##detect NA in Kdis
Kdis.na<-na.omit(Kdis)
miss<-setdiff(Kdis.na[,1],Kdis[,1])
print("diff Kdis na ")
print(miss)
#remove NA in Kdis
for (l in 1:ncol(Kdis)) {ok <- !is.na(Kdis[,l])
                               Kdis[!ok,l]<- mean(Kdis[ok,l])  }
norm<-dcov_mc(Kdis,Kdis,1)
Kdis.list<-as.numeric(c(Kdis))
Kdis.mean<-mean(Kdis.list)
list<- Kdis.list<2*Kdis.mean
print(norm)
print(Kdis.mean)
#for all pheno 
for (i in 1:4) 
{ 
Y<-pheno[,i]
print('i=')
print(i)
Y.produit<-Y %*%t(Y)
Y.produit.list<-as.numeric(c(data.matrix(Y.produit)))

##compute the heritability-linear regression

Kdis.list.red<-as.numeric(as.vector(Kdis.list[list]))
Y.produit.list.red<-as.numeric( as.vector( Y.produit.list[list]))
heri.lm.red<-lm(Y.produit.list.red~Kdis.list.red)
heri.lm.all<-lm(Y.produit.list~Kdis.list)
heri.red<-summary(heri.lm.red)$coef[2,1]
pv.heri.red<-summary(heri.lm.red)$coef[2,4]
heri.all<-summary(heri.lm.all)$coef[2,1]
pv.heri.all<-summary(heri.lm.all)$coef[2,4]
##compute the dcov 
#ny<-length(Y)
#Ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Yp
print(norm)
heri.dcov<-dcov_mc(Y.produit,Kdis,1)
heri.dcov<-heri.dcov/norm
#ne marche pas
#norm.dcov.energy<-(dcov(as.distance(Kdis),as.distance(Kdis),1))^2
#heri.dcov.energy<-(dcov(as.distance(Y.produit),as.distance(Kdis),1))^2
#heri.dcov.energy<-her.dcov.energy/norm.dcov
heri.dcov.energy<-1
norm.dcov<-1
##the results
herita<-c(j,i,heri.red,pv.heri.red,heri.all,pv.heri.all,heri.dcov, heri.dcov.energy, Kdis.mean,norm,norm.dcov)
res<-rbind(res,herita)


}
}
names(res)<-c('pheno', 'method', 'lm.red','pv.lm.red', 'lm','pv.lm','heri.dcov','heri.dcov.energy','meanK', 'dcov(meanK)','dcov(meanK).energy')
rownames(res)<c('K.ibs','K.idb','K.ibs','K.cgta','H.ibs','H.idb','H.ibs','H.cgta','L.ibs','L.idb','L.ibs','L.cgta','T.ibs','T.idb','T.ibs','T.cgta')
write.table(res,file="~/NFG/ler/res.heri.RData", row.names = TRUE,
            col.names = TRUE)









