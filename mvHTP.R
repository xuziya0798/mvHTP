# mvHTP.R
# Hard Thresholding Pursuit algorithm applied in MVMR individual data
# Select the valid IVs and estimate multiple treatments effects. 
#
# Usage: [Shat,Vhat,betahat,ci,betasd] = mvHTP(Y,D,Z,X,s,intercept,alpha,tuning,V,S)
#
# Y: Nx1 outcome vector
# D: Nxq treatment matrix
# Z: Nxpz candidate IVs
# X: Nxpx corvariates
# s: sparsity level
# intercept: whether or not introduce a constant in regression
# alpha: confidence level
# tuning: parameter to adjust threshoulding level in first stage
# oracle: whether or not know the the true S and V
# V: true valid IVs
# S: true relevant IVs
#
# Shat: estimated relevant IVs
# Vhat: estimated valid IVs
# betahat: estimated treatments effects
# ci: alpha-level confidence intervals for treatments effects
# betasd: standard deviation of estimated treatments effects
#


mvHTP <- function(Y,D,Z,X,s,intercept=FALSE,alpha=0.05,tuning=30,oracle=FALSE,V=NULL,S=NULL){
  # Check and Clean Input Type #
  # Check Y
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  
  
  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),(is.matrix(D) || is.data.frame(D)) && ncol(D) >= 1)
  stopifnot(all(!is.na(D)))
  
  
  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),is.matrix(Z))
  stopifnot(all(!is.na(Z)))
  
  # Check dimesions
  stopifnot(length(Y) == nrow(D), length(Y) == nrow(Z))
  
  # Check s
  stopifnot(!missing(s),is.numeric(s) )
  
  # Check X, if present
  if(!missing(X) && !is.null(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))
    
    W = cbind(Z,X)
  } else {
    W = Z
    X = NULL
  }
  
  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)

  # Derive Inputs for TSHT
  n=length(Y)
  pz=ncol(Z)
  p=ncol(W)
  q=length(D)/n
  if(oracle==TRUE){
    stopifnot(!missing(S) && !missing(V))
    Vhat=V
    Shat=S
    
  }else{
    # Estimate Valid IVs
    SetHats = mvHTP.Vhat(Y,D,W,pz,method,intercept,relevant,tuning)
    Vhat = SetHats$Vhat
    Shat = SetHats$Shat
    
  }
  
  # Obtain 2SLS est, se, and ci
  
  W=cbind(Z,X)
  X_=cbind(D,X,Z[,-Vhat])
  auxreg <-lm.fit(W, X_)
  Xhat = as.matrix(auxreg$fitted.values)
  fit <- lm.fit(Xhat,Y)
  ok <- which(!is.na(fit$coefficients))
  yhat <- drop(X_[, ok, drop = FALSE] %*% fit$coefficients[ok])
  res <- Y - yhat
  s2 <- sqrt(sum(res^2)/fit$df.residual)
  betatilde = (fit$coefficients)[1:q]
  
  
  
  H=(s2*solve(t(Xhat)%*%Xhat))[1:q,1:q]
  sd=sqrt(diag(as.matrix(H)))
  ci=cbind(betatilde-qnorm(1-0.05/2)*sd,betatilde+qnorm(1-0.05/2)*sd)
  
  
  return(list(Shat=Shat, Vhat=Vhat, betahat=betatilde, ci=ci, betasd=sd))
  
}

mvHTP.Vhat <- function(Y,D,W,pz,method,intercept=FALSE,relevant,tuning) {
  # Include intercept
  if(intercept) {
    W = cbind(W,1)
  }
  p =  ncol(W) 
  n = length(Y)
  q = length(D)/n
  pj = p+q-1
  
  # ITT effects (OLS Estimation)
  gammatilde=gammahat=matrix(nrow =pz, ncol = q)
  Gammahat=rep(0,pz)
  deltahat=matrix(nrow = n,ncol = q)
  ehat=rep(0,n)
  qrW = qr(W)
  Gammahat= qr.coef(qrW,Y)[1:pz]
  ehat=qr.resid(qrW,Y)
  for(j in 1:q){
    gammahat[,j]=qr.coef(qrW,D[,j])[1:pz]
    deltahat[,j]=qr.resid(qrW,D[,j])
  }
  
  
  # compute the covariance of W,delta(residual of D regress on Z)
  
  Sigmahat=1/(n-p)*t(W)%*%W
  Omegahat22=1/(n-q)*t(deltahat)%*%deltahat
  
  
  
  
  #==========threshold gamma===============
  
    thresh=sqrt(diag(solve(Sigmahat))[1:pz] %*% t(diag(Omegahat22)))*sqrt(tuning*log(n)/n)
    Sflag=(abs(gammahat)>thresh)
    gammatilde=Sflag*gammahat

  
  #========= estimate S* ==============
  Shat=which(apply(Sflag,1,sum)>0) #RESTRICT: ==0
  
  nS=length(Shat)
  
  # Error check
  if(length(Shat) < q){
    warning("VHat Warning: No enough relevant IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being relevant")
    Shat = 1:pz
  }
  
  #=========== estimate V* ===========
  list = HTP(Gammahat[Shat],gammahat[Shat,],1000)
  supp = (list$S)[-(1:q)]-q
  Vhat = (1:pz)[-supp]
  
  # Error check
  if(length(Vhat) < q){
    warning("VHat Warning: No enough valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all relevant IVs being valid")
    Vhat = Shat
  }
  
  return(list(Vhat = Vhat,Shat=Shat))
  
}