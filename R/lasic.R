library(glmnet)
library(ggplot2)
library(magrittr)
library(dplyr)

lasic <- function(fit,X,y,plotit=FALSE) {
  lambda.seq <- fit$lambda
  n <- length(y)
  sigma2 <- estimate_sigma2(fit,X,y)
  results.df <- NULL
  predictions <- predict(fit,X)
  for(i in 1:length(lambda.seq)) {
    lambda <- lambda.seq[i]
    sse <- sum( (predictions[,i] - y)^2 )
    df <- fit$df[i]
    bic.lambda <- BIC(sse,df,sigma2,n)
    aic.lambda <- AIC(sse,df,sigma2,n)
    
    results.df <- rbind(results.df,c(lambda,df,"BIC",bic.lambda))
    results.df <- rbind(results.df,c(lambda,df,"AIC",aic.lambda))
  }
  results.df <- data.frame(results.df,stringsAsFactors=FALSE)
  names(results.df) <- c("Lambda","DF","Type","Value")
  results.df$Lambda <- as.numeric(results.df$Lambda)
  results.df$DF <- as.numeric(results.df$DF)
  results.df$Type <- as.factor(results.df$Type)
  results.df$Value <- as.numeric(results.df$Value)
  best.bic <- results.df %>% filter(Type=="BIC") %>% filter(Value == min(Value))
  best.aic <- results.df %>% filter(Type=="AIC") %>% filter(Value == min(Value))
  if(plotit) {
    lbl <- paste0("BIC Nonzero: ",best.bic$DF,"\nAIC Nonzero: ",best.aic$DF)
    p <- ggplot(results.df,aes(x=log(Lambda),y=Value,color=Type)) + 
          geom_line(size=1) + 
          annotate("text",label=lbl,x=log(median(results.df$Lambda)),y=max(results.df$Value)) +
          theme_bw()
    print(p)
  }
  return(list(AIC=list(Value=best.aic$Value,DF=best.aic$DF,Lambda=best.aic$Lambda),
              BIC=list(Value=best.bic$Value,DF=best.bic$DF,Lambda=best.bic$Lambda)))
}

BIC <- function(sse,df,sigma2,n) {
  return( sse/(n*sigma2) + (log(n)/n)*df )
}

AIC <- function(sse,df,sigma2,n) {
  return( sse/(n*sigma2) + (2/n)*df )
}

estimate_sigma2 <- function(fit,X,y) {
  ## Get sigma2 estimate from the full model ##
  lambda.min <- fit$lambda[which.min(fit$lambda)]
  df.max <- fit$df[which.min(fit$lambda)]
  n <- length(y)
  yhat <- predict(fit,X,s=lambda.min)
  RSS <- sum( (y-yhat)^2 )
  s2 <-  RSS / n ## We're going to use the MLE here
  return( s2 )
}

