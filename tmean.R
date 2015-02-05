meanc=function(X) {
    X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
    return(X)}

dcov_m=function(X,Y) {
    X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
    Y<-t(t(Y)-colMeans(Y))-rowMeans(Y)+mean(Y)
    dcov<- mean(X*Y)
    return(dcov)}
dcor.perso <- function (x, y, index = 1) 
{
  #     if (!(class(x) == "dist")) 
  #         x <- dist(x)
  #     if (!(class(y) == "dist")) 
  #         y <- dist(y)
  #     x <- as.matrix(x)
  #     y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m) 
    stop("Sample sizes must agree")
  if (!(all(is.finite(c(x, y))))) 
    stop("Data contains missing or infinite values")
  if (index < 0 || index > 2) {
    warning("index must be in [0,2), using default index=1")
    index = 1
  }
  stat <- 0
  dims <- c(n, ncol(x), ncol(y))
  Akl <- function(x) {
    d <- as.matrix(x)^index
    rowM <- rowMeans(d)
    colM <- colMeans(d)
    M <- mean(d)
    
    a <- sweep(d, 1, rowM)
    b <- sweep(a, 2, colM)
    
    return(b + M)
  }
  A <- Akl(x)
  B <- Akl(y)
  dCov <- (mean(A * B))
  dVarX <- (mean(A * A))
  dVarY <- (mean(B * B))
  V <- sqrt(dVarX * dVarY)
  if (V > 0) 
    dCor <- dCov/V
  else dCor <- 0
  return(list(dCov = dCov, dCor = dCor, dVarX = dVarX, dVarY = dVarY))
}

covn=function(X,Y){
    return( mean(c(X)*c(Y))) }


resall<-c()
for (p in 1:1000)
{res<-c()
n<-30
#the trait
t<-rnorm(n,mean=0,sd=1)
#the product of the trait
d<-t%*%t(t)
m<-cbind( rbind(d,d), rbind(d,d))

val<-rnorm(n,mean=0,sd=1)
val<-c(val)
relap<-val%*%t(val)
rela<-matrix(rep(val^2,n),ncol=n)+t(matrix(rep(val^2,n),ncol=n))-2*relap
rela<-1-rela
rela.diag<-diag(rela)
rela.off<-rela-rela.diag

#diff<-rnorm(n*n,mean=3,sd=0.2)
#diff<-rep(1,n*n)
#diff<-matrix(diff,ncol=n)
#diag(diff)<-rep(0,n)


for (c in 0:10)
{c<-c*.01
A<-cbind(rela,matrix(rep(0,n*n),ncol=n))
B<-cbind(matrix(rep(0,n*n),ncol=n),rela.diag+rela.off*(1+c))
relat<-rbind(A,B)
res<-cbind( res, covn(m,relat)/sqrt(covn(m,m)*covn(relat,relat)) )

}


for (c in 0:10)
{c<-c*.01
    A<-cbind(rela,matrix(rep(0,n*n),ncol=n))
    B<-cbind(matrix(rep(0,n*n),ncol=n),rela.diag+rela.off*(1+c))
    relat<-rbind(A,B)
    res<-cbind ( res, dcov_m(m,relat)/sqrt(dcov_m(m,m)*dcov_m(relat,relat)) )
    }
res<-cbind(p,res)
resall<-rbind(resall,res)
}
list<-resall[,2]>0
resall.pos<-resall[list,]
mean ((resall.pos[,12]-resall.pos[,2])/resall.pos[,2])
mean ((resall.pos[,7]-resall.pos[,2])/resall.pos[,2])
mean ((resall.pos[,23]-resall.pos[,13])/resall.pos[,13])
mean  ((resall.pos[,18]-resall.pos[,13])/resall.pos[,13])
ex.resall<-resall
resall<-resall.pos

resall.nor<- data.frame(var = resall[,1],
init = resall[,2],
X2= resall[,3],X3= resall[,4],X4= resall[,5],X5= resall[,6],X6= resall[,7],X7= resall[,8],X8= resall[,9],X9=resall[,10],X10= resall[,11], X11=resall[,12],
init.d=resall[,13],
X2.d= resall[,14],X3.d= resall[,15],X4.d= resall[,16],X5.d= resall[,17],X6.d= resall[,18],X7.d= resall[,19],X8.d= resall[,20],X9.d= resall[,21],X10.d= resall[,22],X11.d= resall[,23] )
func<-function(x, y) {x/y}
resall.nor[,2:12]<-apply(resall.nor[,2:12], 2, function(x) func(x,resall.nor[,2])  )
resall.nor[,13:23]<-apply(resall.nor[,13:23], 2, function(x) func(x,resall.nor[,13])  )

res.pet<-resall.nor[1:100,1:12]
res.d.pet<-resall.nor[1:100,c(1,13:23)]

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
resall.melt<-melt(resall.nor, measure.vars=c("init","X2","X3","X4","X5","X6","X7","X8","X9","X10", "X11",
"init.d","X2.d","X3.d","X4.d","X5.d","X6.d","X7.d","X8.d","X9.d","X10.d","X11.d"))

res.sho.melt<-melt(res.pet, measure.vars=c("init","X2","X3","X4", "X5","X6","X7","X8","X9","X10", "X11" ))

res.d.sho.melt<-melt(res.d.pet, measure.vars=c("init.d","X2.d","X3.d","X4.d","X5.d","X6.d","X7.d","X8.d","X9.d","X10.d","X11.d"))

#resall.melt$var <- jitter(resall.melt$var, factor = 0.4)


p <- ggplot(data = res.sho.melt, aes_string(x = "var" , y = "value"))
p<-p+
#geom_line()
geom_point(aes(color = variable, shape =variable),size = 3)+
geom_line(aes(group=var,color=variable),alpha = 0.6, size=0.3)+
#geom_line(aes(linetype = variable))+
#geom_line(aes_string(color ="variable",group="var"),alpha = 0.6, size = 0.8)+
xlab("essai")+
ylab("C effect")+
#scale_shape_manual(values = c(0:1000))+
scale_y_continuous(limits= c(.6,1))+
GetCustomGgplotTheme()
ggsave(plot = p,filename ="effectkin.pdf",width = 11, height = 7)


p <- ggplot(data = res.d.sho.melt, aes_string(x = "var" , y = "value"))
p<-p+
#geom_line()
geom_point(aes(color = variable, shape =variable),size = 3)+
geom_line(aes(group=var,color=variable),alpha = 0.6, size=0.3)+
#geom_line(aes(linetype = variable))+
#geom_line(aes_string(color ="variable",group="var"),alpha = 0.6, size = 0.8)+
xlab("essai")+
ylab("C effect")+
#scale_shape_manual(values = c(0:1000))+
scale_y_continuous(limits= c(.6,1))+
GetCustomGgplotTheme()
ggsave(plot = p,filename ="effectkin.d.pdf",width = 11, height = 7)



