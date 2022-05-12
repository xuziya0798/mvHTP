cv <- function(Y,D,Z,X,lambdaSeq=NA,K = 10,intercept=FALSE) {
  # Include intercept
  if(!missing(X) && !is.null(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))
    
    W = cbind(Z,X)
  } else {
    W = Z
    X = NULL
  }
  if(intercept) {
    W = cbind(W,1)
  }
  
  p =  ncol(W) 
  pz=  ncol(Z)
  n = length(Y)
  q = length(D)/n
  
  if(missing(lambdaSeq) || all(is.na(lambdaSeq))) {
    lambdaSeq=0:(pz-q)
    
  }
  if(any(is.na(lambdaSeq))) {
    warning("Some lambda values are missing. Ignoring these lambda values for cross-validation")
    lambdaSeq = lambdaSeq[!is.na(lambdaSeq)]
  }
  if(length(lambdaSeq) < 2) stop("Only one lambda provided. Please provide multiple lambdas")
  
  lambdaSeq = sort(lambdaSeq,decreasing=TRUE)
  
  
  
  
  # Cross validation
  sampleIndex = sample(rep(1:K,length.out= n))
  errormat = matrix(0,K,length(lambdaSeq))
  for(i in 1:K) {
    testSet = (sampleIndex == i)
    Y.test = Y[testSet]; D.test = D[testSet,,drop=FALSE]; W.test = W[testSet,,drop=FALSE]
    
    testfit = mvHTP.Vhat(Y.test, D.test, W.test,pz=pz,OutputRes=FALSE, intercept = intercept)
     for(j in 1:length(lambdaSeq)){
      trainfit = mvHTP.Vhat(Y[!testSet], D[!testSet,,drop=FALSE], Z[!testSet, , drop=FALSE],OutputRes=FALSE,pz=pz,s=lambdaSeq[j], intercept = intercept)
      residTest = ((testfit$Gammahat) - as.numeric(testfit$gammahat %*% (trainfit$betahat)) - (trainfit$pihat))
      errormat[i,j] = sum(residTest^2)
     }
  }
  cv = colMeans(errormat)
  
  if(all(is.nan(cv))) {
    warning("All lambdas were invalid. Please try different values of lambda for cross-validation to work")
    return(list(lambda = rep(NA,length(lambdaSeq)),estCVError = NA, alpha = rep(NA,ncol(Z)),beta = NA, whichInvalid = NA))
  } else {
    stderror = apply(errormat,2,function(x){sd(x)/sqrt(K)})
    mincv.index = which.min(cv); onestderrorbound = cv[mincv.index] + stderror[mincv.index]
    onestdLambda = max(lambdaSeq[which( (cv <= onestderrorbound) & (cv >= cv[mincv.index]))]) #this is never empty vector 
    
    return(list(s = onestdLambda, estCVError = cv[which(onestdLambda == lambdaSeq)]))
  }
  }
