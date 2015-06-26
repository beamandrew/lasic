library(glmnet)
source('lasic.R')

nobs <- 10000
beta0 <- c(1.5,1.5,0,0,3,0,0,0,0,0)
np <- length(beta0)

X <- matrix(rnorm(nobs*np),nobs,np)

epsilon <- rnorm(nobs,sd=1)
y <- X %*% beta0 + epsilon
y <- c(y)
#X <- scale(X)

fit <- glmnet(x=X,y=y,alpha=1,intercept=FALSE,standardize = FALSE)
params <- lasic(fit,X,y,T)
