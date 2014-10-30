library(plyr)
library(energy)
library(car)
library("matrixStats")
source("mondcov.R")

tfam<-read.table("~/NFG/raw/NFBC_transpose.tfam")
names(tfam)<-c("SUBJID","x","y","t","z","w")

##pheno  in the same order as tfam                           
#testpheno= read.table('~/NFG/raw/pheno2', header=T, sep='')
                

#pheno<-testpheno[,c(-3,-4,-5,-6,-7,-8)]                                       
#opheno<-join(tfam[,c(1,2)],pheno[,c(2,3,4,5,6)],"SUBJID")                       
##for all trait

#Hpheno<-data.matrix(opheno[,c(1,2,4)])                                        
#write.table(Hpheno,file="~/NFG/raw/Hpheno",row.names=F,col.names=F)
#Lpheno<-data.matrix(opheno[,c(1,2,5)])                                        
# write.table(Lpheno,file="~/NFG/raw/Lpheno",row.names=F,col.names=F)
#Tpheno<-data.matrix(opheno[,c(1,2,6)])
# write.table(Lpheno,file="~/NFG/raw/Tpheno",row.names=F,col.names=F)    

# Hpheno<- read.table("~/NFG/raw/Hpheno")
# Lpheno<-read.table("~/NFG/raw/Lpheno")
# Tpheno<-read.table("~/NFG/raw/Tpheno")
# Hpheno<-data.matrix(Hpheno)
# Lpheno<-data.matrix(Lpheno)
# Tpheno<-data.matrix(Tpheno)

##Keep the pheno you wants normalize pheno
#V<-matrix(c(Hpheno[,3],Lpheno[,3],Tpheno[,3]),ncol=3)
#tV<-t(V)-colMeans(V)
#V<-t(tV)
#Sd<-colSds(V)
#tV<-t(V)*(1/Sd)
#V<-t(tV)
#write.table(V,file="~/NFG/raw/phenonor")
V=read.table("~/NFG/raw/phenonor")
source("~/NFG/ler/loadK")
K.divers<-load_K()
typeof(K.divers)
#identical(K_GCTA$id$V1,tfam$SUBJID)


res<-c()
##for all method
for (j in 1:3) 
{Kdis<-K.divers[[j]]
typeof(Kdis)
Kdis<-data.matrix(Kdis)
Ke<-as.numeric(Kdis)
Ke<-c(Ke)
Ke<as.numeric(Ke)

##for all pheno
for (i in 1:3) 
{ Y<-V[,i]
print(i)
Yp<-Y %*%t(Y)
Ye<-data.matrix(Yp)
Ye<-c(Ye)
##compute the heritability-linear regression
Ki<-mean(Ke)
list<- Ke<4*Ki
Kr<-Ke[list]
Yr<-Ye[list]
print(length(Yr))
print(lenght(Ye))
#herlm<-lm(Yr~Kr)
herlm<-0
herlmall<-lm(Ye~Ke)
her<-summary(herlm)$coef[2,1]
pvher<-summary(herlm)$coef[2,4]

##compute the dcov 
ny<-length(Y)
Ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Yp
hercov2<- (0.5)*( ((dcov_m(Ydis2,Kdis,1))^2)/((dcov_m(Kdis,Kdis,1))^2))
hercovp<-  ((dcov_m(Yp,Kdis,1))^2)/ ( (dcov_m(Kdis,Kdis,1))^2)
##the results
phenores<-c(i,j,her,herlmall,herlm,pvher,Ki,hercov2,hercovp)
res<-rbind(res,phenores)
}
}
write.table(res,file="~/NFG/ler/resall.RData")









