pas<-5
for (chr in c(1,4))
{ptdv=read.table(paste("~/NFG/result/Hdap",chr,"pperm",pas,sep=""))


## dessiner les pvalues
for (pheno in 3:6 )
{
jpeg(file=paste("~/NFG/result/Hd",pheno,"-",chr,"pperm",pas,".jpg",sep=""))
plot(ptdv[,2]/1000000,-log(ptdv[,pheno],10) ,col="red",pch=18,xlim=c(40,90),ylim=c(0,13),ann=FALSE)
title(xlab=paste("Chr",chr,"position(Mb) by dc-per pas=",pas,"pheno=",pheno,"Replicate=10^6"), col.lab=rgb(0,0.5,0) )
title(ylab=paste(names(ptdv)[pheno],"-log10(p-vlue)",sep=""),col.lab=rgb(.5,0,0))
dev.off()
}
}

