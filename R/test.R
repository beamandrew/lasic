source('lasic.R')
library(glmnet)
library(dplyr)
library(ggplot2)
## Run tests on lasic and make sure it's working properly ##
p.seq <- c(10,100,1000,5000)
n.reps <- 10
n.seq <- c(100,1000,5000)
frac.nonzero.seq <- c(0.1,0.25,0.5,0.75,1)

results <- NULL
for(p in p.seq) {
  print(p)
  for(frac in frac.nonzero.seq) {
    for(n in n.seq) {
      for(rep in 1:n.reps) {
        X <- matrix(rnorm(sd=1,n=n*p),nrow=n,ncol=p)
        nz <- round(frac*p)
        beta <- rep(0,p)
        beta[1:nz] <- 1
        y <- X %*% beta + rnorm(sd=1,n=n)
        fit <- glmnet(x=X,y=y,alpha=1,intercept=FALSE,standardize = TRUE)
        params <- lasic(fit,X,y,FALSE)
        aic.beta <- coef(fit,s=params$AIC$Lambda)[-1,1] ## Leave out the intercept
        bic.beta <- coef(fit,s=params$BIC$Lambda)[-1,1]
        aic.tp <- length(which( (aic.beta * beta) > 0))/nz
        bic.tp <- length(which( (bic.beta * beta) > 0))/nz
                    
        results <- rbind(results,c(p,frac,nz,n,rep,params$AIC$DF,params$BIC$DF,aic.tp,bic.tp))
      }
    }
  }
}
results <- as.data.frame(results)
names(results) <- c("p","frac_nonzero","num_nonzero","n","rep","AIC_nz","BIC_nz","AIC_tp","BIC_tp")

df <- results %>% 
        group_by(n,p,frac_nonzero,num_nonzero) %>% 
        summarize(AIC_nz=mean(AIC_nz),BIC_nz=mean(BIC_nz),AIC_tp=mean(AIC_tp),BIC_tp=mean(BIC_tp)) %>%
        melt(id=c("n","p","frac_nonzero","num_nonzero"))


df %>% 
  filter(variable == "AIC_tp") %>%
  ggplot(aes(x=n,y=value,color=factor(p))) +
  geom_line() +
  ylab("Fraction of True Covariates with Non-zero Estimates") +
  ggtitle("AIC")

df %>% 
  filter(variable == "BIC_tp") %>%
  ggplot(aes(x=n,y=value,color=factor(p))) +
  geom_line() +
  ylab("Fraction of True Covariates with Non-zero Estimates") +
  ggtitle("BIC")




