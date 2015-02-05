setwd("~/Documents/data/sherlock.nov")
#the lenght of the window for dcov
#pas<-0
pas<-5
#chromosome number
#the phenotype
##the boothstrap parmeter
Repli<-1000000

library(plyr)
library(grid)
datetime.stamp <- format(Sys.time(), "%d%m%Y_%H%M%S")

##load all data
load.chr<-function(chr,pheno){
pv.stand=read.table(paste("resultv/Hdap",chr,"pstand",sep=""),header=T,sep=" ")
pv.stand<-pv.stand[,c(2,pheno)]
names(pv.stand)<-c("posi","pv.stand")  
pv.boot.sel.pas5=read.table(paste("result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",5,"select",sep=""))
pv.boot.sel.pas5<-pv.boot.sel.pas5[,c(3,8)]
names(pv.boot.sel.pas5)<-c("posi","pv.boot.sel.pas5")
if(!chr==4){
pv.boot.sel.pas10=read.table(paste("result/allessaiHdap",chr,"pR","repli=",Repli,"pas=",10,"select",sep=""))
pv.boot.sel.pas10<-pv.boot.sel.pas10[,c(3,8)]
names(pv.boot.sel.pas10)<-c("posi","pv.boot.sel.pas.10")  }
pv.cov.second.0<-read.table(paste("result/allessaiHdap",chr,"pcvsecond0",sep=""), header=T, sep='')
pv.cov.second.0<-pv.cov.second.0[,c(2,pheno)]
names( pv.cov.second.0)<-c("posi","pv.cov.second.0")  
dcor.5=read.table(paste("result/allessaiHdap",chr,"dr",5,sep=""),header=T,sep=" ")
dcor.5<-dcor.5[,c(2,pheno)]
dcor.10=read.table(paste("result/allessaiHdap",chr,"dr",10,sep=""),header=T,sep=" ")
dcor.10<-dcor.10[,c(2,pheno)]
#dcor.15=read.table(paste("result/allessaiHdap",chr,"dr",15,sep=""),header=T,sep=" ")
#dcor.15<-dcor.15[,c(2,pheno)]
pv.dc.ttest.20=read.table(paste("result/allessaiHdap",chr,"pttestpas=",20,sep=""),header=T,sep=" ")
pv.dc.ttest.20<-pv.dc.ttest.20[,c(2,pheno)]
#pv.dcor=read.table(paste("result/Hdap",chr,"pv",10,sep=""),header=T,sep=" "))
#pv.dcor.10<-pv.dcor.10[,c(2,pheno)]
if(chr==19){
pv.dc.ttest.20.full=read.table(paste("result/allessaiHdap",chr,"fullpttestpas=",20,sep=""),header=T,sep=" ")
pv.dc.ttest.20.full<-pv.dc.ttest.20.full[,c(2,pheno)]
pv.dc.ttest.10.full=read.table(paste("result/allessaiHdap",chr,"fullpttestpas=",10,sep=""),header=T,sep=" ")
pv.dc.ttest.10.full<-pv.dc.ttest.10.full[,c(2,pheno)]
pv.dc.ttest.5.full=read.table(paste("result/allessaiHdap",chr,"fullpttestpas=",5,sep=""),header=T,sep=" ")
pv.dc.ttest.5.full<-pv.dc.ttest.5.full[,c(2,pheno)] } 
#pv.dc.moi.10=read.table(paste("result/allessaiHdap",chr,"pdrmoi",10,sep=""),header=T,sep=" ")
#pv.dc.moi.10<-pv.dc.moi.10[,c(2,pheno)]
pv.dc.moi.20=read.table(paste("result/allessaiHdap",chr,"pdrmoi",20,sep=""),header=T,sep=" ")
pv.dc.moi.20<-pv.dc.moi.20[,c(2,pheno)]
dc.moi.10=read.table(paste("result/allessaiHdap",chr,"drmoi",10,sep=""),header=T,sep=" ")
dc.moi.10<-dc.moi.10[,c(2,pheno)]
#creation of data
res<-pv.boot.sel.pas5 
names(res)<-c("posi","pv.boot.sel.pas5")  
if (!chr==4){
  res<-join(pv.boot.sel.pas5,pv.boot.sel.pas10,by="posi",type="full")
  names(res)<-c("posi","pv.boot.sel.pas5","pv.boot.sel.pas10")  
  }
#to order p.stand as the other ones
res.int<-join(pv.cov.second.0,pv.stand,by="posi",type="full")
#remove pv.cov.second.0
res.int<-res.int[,c(1,2)]
names(res.int)<-c('posi','pv.stand')
res.fin<-cbind( #pv.cov.second.0
  #dcor.5[,-1],dcor.10[,-1],
  #pv.dcor.10[,-1],
  pv.dc.ttest.20
  #dc.moi.10[,-1],
  #pv.dc.moi.10[,-1],
  #pv.dc.moi.20[,-1] 
  )
res.int<-join(res.int,res.fin,by="posi",type="full") 
  debut.names<-c('posi','pv.stand',
                  #'pv.cov.sec',
                  #'dcor.5','dcor.10',
                  #'pv.dcor.10',
                  'pv.dc.ttest.20'
                  #'dc.moi.10',
                  #'pv.dc.moi.10,
                  #'pv.dc.moi.20'
                 )
names(res.int)<-debut.names
if(chr==19){
  toto<-cbind(pv.dc.ttest.5.full ,pv.dc.ttest.10.full[,-1] #,pv.dc.ttest.20.full[,-1]
              )
    res.int<-join(res.int, toto,by="posi",type="full") 
    names(res.int)<-c(debut.names, 'pv.dc.ttest.5.full'
                      ,'pv.dc.ttest.10.full'
                      #,
                      #'pv.dc.ttest.20.full'
                      )
    }
res<-join(res,res.int,by="posi",type="full")
##format for picture
fonc<-function(x) {x/1000000}
res[,1]<-sapply(res[,1], function(x) fonc(x)  )
func<-function(x) {-log(x,10)}
res[,2:ncol(res)]<-apply(res[,2:ncol(res)], 2, function(x) func(x)  )

##selection draw
debut.var<- c('pv.dc.ttest.20','pv.stand', 'pv.boot.sel.pas5'
                                 #'pv.cov.sec',
                                   #'dcor.5','dcor.10',
                                 #'pv.dc.ttest.5',
                                   #'dc.moi.10',
                                   #'pv.dc.moi.20'
                                   )
if (!chr==4) debut.var<- c('pv.dc.ttest.20','pv.stand','pv.boot.sel.pas5','pv.boot.sel.pas10'
                          #'pv.cov.sec',
                          #'dcor.5','dcor.10',
                          #'pv.dc.ttest.5',
                          #'dc.moi.10',
                          #'pv.dc.moi.20'
                          )
all.var<-debut.var
if(chr==19) all.var<-c('pv.dc.ttest.5.full',debut.var
                       #'pv.dc.ttest.20.full',
                       'pv.dc.ttest.10.full'
                       )
                      

res.melt<-melt(res, measure.vars= all.var)
return(res.melt)
}


library(ggplot2)
library(reshape2)

GetCustomGgplotTheme <- function() {
  #theme_bw() +
  theme(title = element_text(size = rel(1)),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22, face = "bold"),
        plot.margin = grid::unit(c(1,1,0,0), "cm"),
        legend.justification = c(1,0),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        legend.key.width = grid::unit(0.5, "line"),
        legend.key.height = grid::unit(1.8,"line"),
        legend.key.size = grid::unit(3, "line"))
}
create_ggplot <- function(data,chr) {  pv.stand=read.table(paste("resultv/Hdap",1,"pstand",sep=""),header=T,sep=" ")
                                       pv.stand<-pv.stand[,c(2,pheno)]
  ggplot(data = data, aes_string(x = "posi" , y = "value"))+
  geom_point(aes(color = variable,size=variable),shape=1)+
  xlab("position")+
  ylab(paste( "pv divers", names(pv.stand)[2] )) +
  scale_shape_discrete(solid=T, legend=F) +
  scale_y_continuous(limits= c( 0,20))+
  scale_x_continuous(limits= c( 40,90))+ 
  GetCustomGgplotTheme()+
  ggtitle(paste0("pv value versus ttest.dcov chr=",chr,"pheno=",names(pv.stand)[2]))
}

library(gridBase)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pheno<-4


p<-create_ggplot( load.chr(1,pheno),1 )

q<-create_ggplot(load.chr(4,pheno),4)

r<-create_ggplot(load.chr(19,pheno),19)

pv.stand=read.table(paste("resultv/Hdap",1,"pstand",sep=""),header=T,sep=" ")
pv.stand<-pv.stand[,c(2,pheno)]

pdf(file = paste0("pvalue", "pheno=",names(pv.stand)[2],"-", datetime.stamp, ".pdf"), 
    width = 30, 
    height = 5)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))   

print(p, vp = vplayout(1,1))
print(q, vp = vplayout(1,2))        
print(r, vp = vplayout(1,3))

upViewport(0)

dev.off() 



# to select
#threshold <- 0.15
#df <- data.frame(phenotype.similarity = Y.square.vec[which(phi >= threshold)],
                 #genotype.similarity = phi[which(phi >= threshold)])







 


