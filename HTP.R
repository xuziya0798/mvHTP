# HTP.m
# Basic implementation of the Hard Thresholding Pursuit algorithm
# Find the s-sparse solution of the underdetermined mXN linear system Ax=y 
#
# Usage: [x,S,NormRes,NbIter] = HTP(y,A,s,MaxNbIter,mu,x0,TolRes,Warnings,Eps)
#
# y: mx1 measurement vector
# A: mxN measurement matrix
# s: sparsity level
# MaxNbIter: number of iterations not to be exceeded (optional, default=500)
# mu: step factor of A%*%(y-Ax) in the support-identification step (optional, default=1)
#     use 'NHTP' for the Normalized Hard Thresholding Pursuit algorithm 
# x0: initial vector (optional, default=zero)
# TolRes: tolerance on the Euclidean norm |y-Ax|/|y| of the relative residual (optional, default=1e-4)
# Warnings: use 'On' to display the warnings, 'No' otherwise (optional, default='On')
# Eps: thresholding for sending to zero the entries of a vector with magnitude smaller than Eps times the largest entry of the vector in magnitude (optional, default=1e-8)
#
# x: an s-sparse vector which is the possible solution of Ax=y
# S: the support of x
# NormRes: the Euclidean norm |y-Ax| of the residual
# NbIter: the number of iterations performed until stationary outputs are reached
#         if stationary outputs are not reached, a warning is displayed and NbIter takes MaxNbIter 
#
# Written by Simon Foucart in August 2010, updated in February 2011
# Code proposed and used in the paper "Hard Thresholding Pursuit: an algorithm for Compressive Sensing"
# Send comments to simon.foucart@centraliens.net


HTP <- function(y,A,s,MaxNbIter=500,mu=1,x0,TolRes=1e-4,Warnings='On',Eps=1e-8){

N = ncol(A)
## set the default values
if(x0=NULL){
  x0=rep(0,N)
  S0=1:s
}




## define auxiliary quantities
B = t(A) %*% A 
z = t(A) %*% y 
S0 = order(abs(x0),decreasing=TRUE)[1:s]


## initialization
x=x0 
S=S0  
g=z-B%*%x 
if(mu=='NHTP') {
  Mu=(norm(g[S])/norm(A[,S]%*%g[S]))^2 
}
else  Mu=mu 

v=x+Mu%*%g  
absv=abs(v)  
zero_idx=which(absv<Eps%*%max(absv)) 
absv(zero_idx)=rep(0,length(zero_idx)) 
Snew=order(abs(absv),decreasing=TRUE)[1:s] 
xnew=rep(0,N) 
xnew[Snew]=solve(A[,Snew],y )
NbIter=1 


## main loop
while ( (sum(S==Snew)<s) && (NbIter<MaxNbIter) ){
  x=xnew 
  S=Snew 
  g=z-B%*%x 
  
  if (mu=='NHTP') Mu=(norm(g[S])/norm(A[,S]%*%g[S]))^2 
  else Mu=mu 
 
  v=x+Mu%*%g  
  absv=abs(v)  
  zero_idx=which(absv<Eps%*%max(absv)) 
  absv[zero_idx]=rep(0,length(zero_idx)) 
  
  Snew=order(absv,decreasing=TRUE)[1:s]  
  xnew=rep(0,N) 
  xnew[Snew]=solve(A[,Snew],y )
  NbIter=NbIter+1 
}
  

## outputs
NormRes=norm(y-A%*%xnew)
if(Warnings=='On'){
  if(sum(S==Snew)<s) {
    cat('Warning: HTP did not converge when using a number of iterations =', MaxNbIter,"\n")
  }
  else {
    if (NormRes>TolRes%*%sqrt(sum(y^2))){
      disp(cat('Warning: HTP converged to an incorrect solution (norm of residual =', NormRes,')')) 
    }
  }
}


x=xnew 
S=Snew 
return(list(x=x,
            S=S,
            NormRes=NormRes,
            NbIter=NbIter))

}