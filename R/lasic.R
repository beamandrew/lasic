library(glmnet)

lasic <- function(fit,X,y,metric="BIC",df.method="naive") {
  lambda.seq <- fit$lambda
  n <- length(y)
  sigma2 <- 1
  bic.seq <- NULL
  aic.seq <- NULL
  for(i in 1:length(lambda.seq)) {
    lambda <- lambda.seq[i]
    sse <- sum( (predict(fit,X,s=lambda) - y)^2 )
    df <- fit$df[i]
    bic.lambda <- BIC(sse,df,sigma2,n)
    bic.seq <- c(bic.seq,bic.lambda)
    
    aic.lambda <- BIC(sse,df,sigma2,n)
    aic.seq <- c(aic.seq,aic.lambda)
  }  
}

BIC <- function(sse,df,sigma2,n) {
  return( sse/(n*sigma2) + (log(n)/n)*df )
}

AIC <- function(fit,df) {
  return( sse/(n*sigma2) + (log(n)/n)*df )
}
