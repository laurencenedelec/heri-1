##function to compute the  dcov  matrices form
dcov_m=function(xdis,ydis,index) {
X<-xdis^index
Y<-ydis^index
X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
Y<-t(t(Y)-colMeans(Y))-rowMeans(Y)+mean(Y)
dcov2<- mean(X*Y)
if (dcov2<0) {dcov<- -sqrt(-dcov2)} else {dcov<-sqrt(dcov2)}
return(dcov)}

dcov_mc=function(xdis,ydis,index) {
X<-xdis^index
Y<-ydis^index
#Y<-sqrt(Y)
X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
Y<-t(t(Y)-colMeans(Y))-rowMeans(Y)+mean(Y)
dcov<- mean(X*Y)
return(dcov)}


##function to compute the  dcov  vector matrix mix form
dcov_vm=function(x,ydis,index) {
X<-data.matrix(x)
nx<-length(X)
Xe<-data.matrix(X%*%t(X))
xdis2<- matrix(rep(X^2,nx),ncol=nx)+t(matrix(rep(X^2,nx),ncol=nx))-2*Xe
xdis<-xdis2^(1/2)
ydis<-data.matrix(ydis)
dcov<-dcov_m(xdis,ydis,index)
return(dcov)}


##function to compute the  dcor  vector matrix mix form
dcor_vm=function(x,ydis,index) {
X<-data.matrix(x)
nx<-length(X)
Xe<-data.matrix(X%*%t(X))
xdis2<- matrix(rep(X^2,nx),ncol=nx)+t(matrix(rep(X^2,nx),ncol=nx))-2*Xe
xdis<-xdis2^(1/2)
dcor<-dcov_m(xdis,ydis,index)/sqrt( dcov_m(xdis,xdis,1)*dcov_m(ydis,ydis,1)  )
return(dcor)}


##function to compute the dcov agit sur les vectors
dcov_vv=function(x,y,index)
{X<-data.matrix(x)
nx<-length(X)
Xe<-data.matrix(X%*%t(X))
xdis2<- matrix(rep(X^2,nx),ncol=nx)+t(matrix(rep(X^2,nx),ncol=nx))-2*Xe
xdis<-xdis2^(1/2)
Y<-data.matrix(y)
ny<-length(Y)
Ye<-data.matrix(Y%*%t(Y))
ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Ye
ydis<-ydis2^(1/2)
dcov<- dcov_m(xdis,ydis,index)
return(dcov)}

##function to compute the dcor agit sur les vectors
dcorper_vv=function(x,y,index)
{X<-data.matrix(x)
nx<-length(X)
Xe<-data.matrix(X%*%t(X))
xdis2<- matrix(rep(X^2,nx),ncol=nx)+t(matrix(rep(X^2,nx),ncol=nx))-2*Xe
xdis<-xdis2^(1/2)
Y<-data.matrix(y)
ny<-length(Y)
Ye<-data.matrix(Y%*%t(Y))
ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Ye
ydis<-ydis2^(1/2)
dcor<- dcov_m(xdis,ydis,index)/(sqrt(dcov_m(xdis,xdis,index)*dcov_m(ydis,ydis,index)))
return(dcor)}

##function to compute a test: matrix from
pcov_m=function(xdis,ydis,alpha)
{ptest<- ncol(xdis)*( (dcov_m(xdis,ydis,1))^2/(mean(xdis)*mean(ydis)))
return(sqrt(ptest)-qnorm(1-alpha/2))}

##function to compute a test : vector form dimension one
pcov_v=function(x,y)
{X<-data.matrix(x)
nx<-length(X)
Xe<-data.matrix(X%*%t(X))
xdis2<- matrix(rep(X^2,nx),ncol=nx)+t(matrix(rep(X^2,nx),ncol=nx))-2*Xe
xdis<-xdis2^(1/2)
Y<-data.matrix(y)
ny<-length(Y)
Ye<-data.matrix(Y%*%t(Y))
ydis2<- matrix(rep(Y^2,ny),ncol=ny)+t(matrix(rep(Y^2,ny),ncol=ny))-2*Ye
ydis<-ydis2^(1/2)
ptest<- ncol(xdis)*((dcov_m(xdis,ydis,1))^2/(mean(xdis)*mean(ydis)))
return(2-2*pnorm(sqrt(ptest)))
}

##function to compute a test : vector from dimension bigger
pcov_vvdim=function(x,y)
{x<-data.matrix(x)
nx<-nrow(x)
Xe<-data.matrix(x%*%t(x))
xdis2<-matrix(rep(rowSums(x^2),nx),nx)+t( matrix(rep(rowSums(x^2),nx),nx))-2*Xe
xdis<-xdis2^(1/2)
y<-data.matrix(y)
ny<-nrow(y)
Ye<-data.matrix(y%*%t(y))
ydis2<-matrix(rep(rowSums(y^2),ny),ny)+t(matrix(rep(rowSums(y^2),ny),ny))-2*Ye
ydis<-ydis2^(1/2)
ok <- !is.na(ydis)
ydis[!ok]<-0
ptest<- ncol(xdis)*((dcov_m(xdis,ydis,1))^2/(mean(xdis)*mean(ydis)))
if (ptest>0)  {ptest<- 2-2*pnorm(sqrt(ptest) )} else {ptest<-0}
return(ptest)}

##function to compute dcor mix vector and  v dimension bigger
dcor_vvdim=function(x,y)
{x<-data.matrix(x)
nx<-length(x)
Xe<-data.matrix(x%*%t(x))
xdis2<-matrix(rep(rowSums(x^2),nx),nx)+t( matrix(rep(rowSums(x^2),nx),nx))-2*Xe
xdis<-xdis2^(1/2)
y<-data.matrix(y)
ny<-nrow(y)
Ye<-data.matrix(y%*%t(y))
ydis2<-matrix(rep(rowSums(y^2),ny),ny)+t( matrix(rep(rowSums(y^2),ny),ny))-2*Ye
ydis<-ydis2^(1/2)
ok <- !is.na(ydis)
ydis[!ok]<-0
return( dcov_m(xdis,ydis,1)/(sqrt( dcov_m(xdis,xdis,1)*dcov_m(ydis,ydis,1)) ) )
}
