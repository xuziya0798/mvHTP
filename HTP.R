# HTP.R
# Hard Thresholding Pursuit algorithm applied in MVMR summary statistics
# Find the s-sparse solution of the underdetermined px(p+q) linear system Ax=y 
#
# Usage: [x,S,NormRes,NbIter] = HTP(Gamma,gamma,s,MaxNbIter,mu,x0,TolRes,Warnings,Eps)
#
# Gamma: px1 measurement vector
# gamma: pxq measurement matrix
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
# Translate the original matlab code, written by Simon Foucart in August 2010, to R




HTP <- function(Gamma,gamma,s,MaxNbIter=1000,mu=1/5,x0=NULL,TolRes=1e-4,Warnings='On',Eps=1e-8){

p = nrow(gamma)
q = ncol(gamma)
## set the default values
if(is.null(x0)){
  x0=c(rep(1,q),rep(0,p))
  S0=c(1:q, q + 1:s)
}

## renormalization of gamma
d=rep(1,q);
for(j in 1:q){
  d[j]=1/sqrt(sum(gamma[,j]^2))
} 

gamma0=gamma %*% diag(d)

A=cbind(gamma0,diag(1,p))


## define auxiliary quantities
B = t(A) %*% A 
z = t(A) %*% Gamma 
S0 = c(1:q, q+ order(abs(x0[-(1:q)]),decreasing=TRUE)[1:s])


## initialization
x=x0 
S=S0  
g=z-B%*%x 
if(mu=='NHTP') {
  Mu=(norm(g[S])/norm(A[,S]*g[S]))^2 
}
else  Mu=mu 

v=x+Mu*g  
absv=abs(v[-(1:q)]) 
zero_idx=q+which(absv<Eps*max(absv)) 
absv[zero_idx]=rep(0,length(zero_idx)) 
Snew=c(1:q, q + order(abs(absv),decreasing=TRUE)[1:s]) 
xnew=rep(0,p+q) 
xnew[Snew]=solve(t(A[,Snew])%*%A[,Snew],t(A[,Snew])%*%Gamma )
NbIter=1 


## main loop
while ( (sum(S==Snew)<s+q) && (NbIter<MaxNbIter) ){
  x=xnew 
  S=Snew 
  g=z-B%*%x 
  
  if (mu=='NHTP') Mu=(norm(g[S])/norm(A[,S]*g[S]))^2 
  else Mu=mu 
 
  v=x+Mu*g  
  absv=abs(v[-(1:q)]) 
  zero_idx=q+which(absv<Eps*max(absv)) 
  absv[zero_idx]=rep(0,length(zero_idx))
  
  Snew=c(1:q, q + order(abs(absv),decreasing=TRUE)[1:s]) 
  xnew=rep(0,p+q) 
  xnew[Snew]=solve(t(A[,Snew])%*%A[,Snew],t(A[,Snew])%*%Gamma )
  NbIter=NbIter+1 
}
  

## outputs
NormRes=norm(Gamma-A%*%xnew)
if(Warnings=='On'){
  if(sum(S==Snew)<s+q) {
    cat('Warning: HTP did not converge when using a number of iterations =', MaxNbIter,"\n")
  }
  else {
    if (NormRes>TolRes*sqrt(sum(Gamma^2))){
      cat('norm of residual =', NormRes)
    }
  }
}

x=xnew
x[1:q]=diag(d) %*% xnew[1:q]
S=Snew 
return(list(x=x,
            S=S,
            NormRes=NormRes,
            NbIter=NbIter))
}