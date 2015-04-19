
options(digits=16)
#Amat=matrix(c(1,0,1,0,1,0,0,1,1,0,1,0,1,1,0),byrow=T,ncol=3)
Amat=matrix(c(3,1,2,0,3,0,1,2,4,1,2,2,0,4,3),byrow=T,ncol=3)
#Amat=matrix(runif(10,0,1),byrow=T,ncol=2)

Amat

A<-Amat
dimension<-dim(A)
dimension
m<-dimension[1]
n<-dimension[2]
V<-diag(rep(1,n))
V
no_rot=0 #contador de no rotaciones
ban=F
col_ort=n*(n-1)/2 #número de col ortogonales
tol=1.e-8
#tol=.Machine$double.xmin

signo<-function(x){
  ifelse(x>=0,1,-1)
}
sweep=0;
maxsweep=10;

while(!ban && sweep<=maxsweep){
no_rot=0
for(i in 1:(n-1)){
for(j in (i+1):n){
  ai=A[,i]
  aj=A[,j]
  prod=sum(ai*aj)
  ainorma=sqrt(sum(ai*ai))
  ajnorma=sqrt(sum(aj*aj))
  #if(prod==0){
   # costeta=1
    #senteta=0
  #}
  #else{
  num=ajnorma^2-ainorma^2
  den=2*prod
  if(abs(den)<=1.e-16*abs(num)){
    costeta=1
    senteta=0
    no_rot=no_rot+1
  }
  else{
    ji=num/den
    
    t=signo(ji)/(abs(ji)+sqrt(1+ji^2))
    
    costeta=1/sqrt(1+t^2)
    senteta=costeta*t
    cte=prod/(ainorma*ajnorma)
    if(abs(cte)>tol){
      vi=V[,i]
      vj=V[,j]
      mataux=rbind(cbind(ai,aj),cbind(vi,vj))
      rot=matrix(c(costeta,senteta,-senteta,costeta),byrow=T,ncol=2)
      mataux=mataux%*%rot
      #print(mataux)
      A[,i]=mataux[1:m,1]
      A[,j]=mataux[1:m,2]
      V[,i]=mataux[(m+1):(m+n),1]
      V[,j]=mataux[(m+1):(m+n),2]
    }
    else{
      no_rot=no_rot+1
    }
    
  #  }
    }#fin else
  print(no_rot)
  
}#fin for

}#fin for
  sweep=sweep+1

if(no_rot==col_ort){
  ban=T
}

}#fin while
no_rot

sweep

#valores singulares de A, nombres de acuerdo a columnas de A
sigma1=sqrt(sum(A[,1]*A[,1]))
sigma2=sqrt(sum(A[,2]*A[,2]))
sigma3=sqrt(sum(A[,3]*A[,3]))
Sigma=diag(c(sigma1,sigma2,sigma3))
#Sigma=diag(c(sigma1,sigma2))

U=scale(A,center=F,scale=c(sigma1,sigma2,sigma3))
#U=scale(A,center=F,scale=c(sigma1,sigma2))

Amat

U%*%Sigma%*%t(V)

svd(Amat)
U
V
A