m<-c(1,2,0,0,6,8,0,0,0,0,1,2,0,0,6,8)

#m<-c(1,2,4,0,0,1,10,3,0,0,3,5,6,0,0,0,0,0,3,5,0,0,0,13,6)
m<-matrix(m,ncol=4)

#n<-c(1,2,4,0,0,1,10,3,0,0,3,5,6,0,0,0,0,0,3+c,5+c,0,0,0,13+c,6+c)
#n<-matrix(n,ncol=5)
meanc=function(X) {
X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
return(X)}


dcov_m=function(X,Y) {
X<-t((t(X)-colMeans(X)))-rowMeans(X)+mean(X)
Y<-t(t(Y)-colMeans(Y))-rowMeans(Y)+mean(Y)
dcov2<- mean(X*Y)
if (dcov2<0) {dcov<- -sqrt(-dcov2)} else {dcov<-sqrt(dcov2)}
return(dcov)}

covn=function(X,Y){
return(sqrt( mean(c(X)*c(Y))))}

#meanc(X)
for (c in 1:40)
{c<-c*10
 #n<-c(1,2,4,0,0,1,10,3,0,0,3,5,6,0,0,0,0,0,3+c,5+c,0,0,0,13+c,6+c)
 n<-c(3,5,0,0,13,6,0,0,0,0,3+c,5+c,0,0,13+c,6+c)
 n<-matrix(n,ncol=4)
print(dcov_m(m,n)/sqrt(dcov_m(m,m)*dcov_m(n,n)))
print(covn(m,n)/sqrt(covn(m,m)*covn(n,n)) )}
