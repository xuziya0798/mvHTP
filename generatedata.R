generatedata <- function(n,r,beta,gamma,pi,phi=NULL,psi=NULL){
  #=====generate data=============
  set.seed(r)
  q=length(beta)
  pz=length(pi)
  px=length(psi)
  p=px+pz
  
  #Sigma = matrix(nrow =p , ncol = p)
  # for(i in 1:p){
  #   for (j in 1:i) {
  #     Sigma[i,j]=0.5^(i-j)
  #     Sigma[j,i]=Sigma[i,j]
  #   }
  # }
  
  Sigma=diag(p)
  
  library(MASS)
  library(mvtnorm)
  mu=rep(0,p)
  W=mvrnorm(n,mu,Sigma)
  if(px!=0){
    X=W[,pz+1:px]
  } else X=NULL
  Z=W[,1:pz]
  
  # Theta22= matrix(nrow = q,ncol = q)
  # for (i in 1:q) {
  #   for (j in 1:i) {
  #     Theta22[i,j]=0.5^(i-j)
  #     Theta22[j,i]=Theta22[i,j]
  #   }
  # }
  # Theta12=rep(0.5,q)
  # 
  # Theta=cbind(1,t(Theta12))
  # Theta=rbind(Theta,cbind(Theta12,Theta22))
  
  Theta=matrix(0.25,1+q,1+q)+diag(1+q)*0.55
  #Theta=diag(1+q)  ##id case
  mu_e=rep(0,1+q)
  e=mvrnorm(n,mu_e,Theta)
 
  
  delta=e[,1+1:q]
  if(px!=0){
    D=X%*%phi+Z%*%gamma+delta
    Y=D%*%beta+X%*%psi+Z%*%pi+e[,1]
  } else{
    D=Z%*%gamma+delta
    Y=D%*%beta+Z%*%pi+e[,1]
  }
  S=1:pz
  
    S=S[(apply(gamma!=0,1,sum))>0]
    V=S[pi[S]==0]

  
  return(list(Y=Y,Z=Z,X=X,D=D,S=S,V=V))
}
  
  