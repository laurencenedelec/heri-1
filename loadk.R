load_K <- function() {

tfam<-read.table("~/NFG/raw/NFBC_transpose.tfam")
names(tfam)<-c("SUBJID","x","y","t","z","w")

## pourcentage of allele share between two indiv IBS  square matrix
K.ibs= read.table('~/NFG/result/NFGsimi.mibs', sep='')
K.ibs.id=read.table('~/NFG/result/NFGsimi.mibs.id', sep='')
colnames(K.ibs) <- paste(K.ibs.id$V1)
rownames(K.ibs) <- paste(K.ibs.id$V1)
K.ibs<-1-K.ibs
write.table(K.ibs[1:10,1:10],file="~/NFG/raw/K.ibs", row.names = TRUE,
            col.names = TRUE)
print('K.ibs')




##Other K ibd compute from IBS there is one more individu ??
pt=read.table("~/NFG/result/NFGredibsibd.genome" ,header=TRUE, sep='')
K.ibsg<- matrix( 0, nrow = nrow(tfam)+1,ncol = nrow(tfam)+1)
K.ibdg<- matrix(data = 0, nrow = nrow(tfam)+1,ncol = nrow(tfam)+1)
for (j in 1:nrow(pt))
{K.ibsg[pt$IID1[j],pt$IID2[j]]<-pt$PI_HAT[j]
K.ibdg[pt$IID1[j],pt$IID2[j]]<-pt$DST[j]
K.ibsg[pt$IID2[j],pt$IID1[j]]<-pt$PI_HAT[j]
K.ibdg[pt$IID2[j],pt$IID1[j]]<-pt$DST[j]
}
K.ibsg<-K.ibsg[as.vector(tfam$SUBJID),as.vector(tfam$SUBJID)]
K.ibdg<-K.ibdg[as.vector(tfam$SUBJID),as.vector(tfam$SUBJID)]
colnames(K.ibsg) <- paste(tfam$SUBJID)
rownames(K.ibsg) <- paste(tfam$SUBJID)
colnames(K.ibdg) <- paste(tfam$SUBJID)
rownames(K.ibdg) <- paste(tfam$SUBJID)
for (j in 1:(nrow(tfam)))
{K.ibsg[j,j]<-0
K.ibdg[j,j]<-1}
#K.ibsg<-(1-K.ibsg)
K.ibdg<-(1-K.ibdg)

write.table(K.ibsg[1:10,1:10],file="~/NFG/raw/K.ibsg",row.names = TRUE,
            col.names = TRUE)
write.table(K.ibdg[1:10,1:10],file="~/NFG/raw/K.ibdg",row.names = TRUE,
            col.names = TRUE)
#kdisibs<-read.table("~/NFG/raw/K.ibsg")
#kdisibd<-read.table("~/NFG/raw/K.ibdg")

print('K ibs ibd')

##same as K.ibs from plink
#K.ibs.plink<-read.table("~/NFG/result/plinkdis.mdist")
#colnames(K.ibs.plink) <- paste(tfam$SUBJID)
#rownames(K.ibs.plink) <- paste(tfam$SUBJID)




##GCTA K normalize
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}


 # Data are stored in the data/METHODNAME directory
 K_GCTA <- ReadGRMBin(prefix = paste0("~/NFG/raw/Kcgta"))

  # Create the matrix Ki_GCTA from the data
  # Data come as a lower triangular matrix
  Ki_GCTA<- diag(nrow(K_GCTA$id))
  diag(Ki_GCTA) <-1- K_GCTA$diag
  Ki_GCTA[lower.tri(Ki_GCTA, diag = F)] <- K_GCTA$off
  Ki_GCTA <- Ki_GCTA + t(Ki_GCTA) - diag(diag(Ki_GCTA))
  # Give the correct col/row names (this information is used later on to filter data)
  colnames(Ki_GCTA) <- paste(K_GCTA$id$V1)
  rownames(Ki_GCTA) <- paste(K_GCTA$id$V1)

write.table(Ki_GCTA[1:10,1:10],file="~/NFG/raw/kgcta",row.names = TRUE,
            col.names = TRUE)

#identical(K_GCTA$id$V1,tfam$SUBJID)
print('Kgcta')

##GCTA avec plink 1.9
id<- read.table('~/NFG/result/plinkdisex.mibs.id')
Ki_pGC<- diag(length(id$V1))
Ki_pGC_tri<- data.matrix( read.table('~/NFG/result/plinkdisex.mibs', fill=TRUE, col.names=paste("V", 1:length(id$V1))))
Ki_pGC[lower.tri(Ki_pGC, diag = T)] <- Ki_pGC_tri[lower.tri(Ki_pGC_tri,diag=T)]
Ki_pGC <- Ki_pGC + t(Ki_pGC) - diag(diag(Ki_pGC))
rownames(Ki_pGC) <- paste(id$V1)
colnames(Ki_pGC) <- paste(id$V1)
Ki_pGC<-data.matrix(Ki_pGC)
Ki_pGC<-1-Ki_pGC
write.table(Ki_pGC[1:10,1:10],file="~/NFG/raw/kpgc",row.names = TRUE,
            col.names = TRUE)
print('KpGC')

##GCTA with plink -make-rel 
K_ppGC.id<- read.table('~/NFG/result/plinkgcta.rel.id')
K_ppGC<- data.matrix( read.table('~/NFG/result/plinkgcta.rel'))
rownames(K_ppGC) <- paste(K_ppGC.id$V1)
colnames(K_ppGC) <- paste(K_ppGC.id$V1)
write.table(K_ppGC[1:10,1:10],file="~/NFG/raw/kppgc", row.names = TRUE,
            col.names = TRUE)
diag(K_ppGC)<-1-diag(K_ppGC)
print('KppGC')

K.divers<-list(K_ibs=K.ibsg,K_ibd=K.ibdg,K_plink=K.ibs,K_gcta=Ki_pGC,K_ppGC=K_ppGC,K_GCTA=Ki_GCTA)
return(K.divers)

}