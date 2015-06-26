library(glmnet)
source('lasic.R')

nobs <- 1000
np <- 10

X <- matrix(rnorm(nobs*np),nobs,np)
beta0 <- c(1,1.5,0,0,3,0,0,0,0,0)
epsilon <- rnorm(nobs,sd=3)
y <- X %*% beta0 + epsilon
y <- c(y)
X <- scale(X)

fit <- glmnet(x=X,y=y,alpha=1,intercept=FALSE,standardize = TRUE)
lasic(fit,X,y,T)
