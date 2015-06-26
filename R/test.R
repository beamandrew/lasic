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
        summarize(AIC_nz=mean(AIC_nz),
                  BIC_nz=mean(BIC_nz),
                  AIC_tp=mean(AIC_tp),
                  BIC_tp=mean(BIC_tp)) %>%
        mutate(AIC_eff=num_nonzero/AIC_nz,BIC_eff=num_nonzero/BIC_nz) %>%
        melt(id=c("n","p","frac_nonzero","num_nonzero"))


df %>% 
  filter(variable == "AIC_tp" | variable == "BIC_tp") %>%
  ggplot(aes(x=n,y=value,color=variable)) +
  facet_wrap(~ p) +
  geom_line() +
  ylab("Fraction of True Covariates with Non-zero Estimates") +
  theme_bw()

df %>% 
  filter(variable == "AIC_eff" | variable == "BIC_eff") %>%
  ggplot(aes(x=n,y=value,color=variable)) +
  facet_wrap(~ p) +
  geom_line() +
  ylab("Effeciency") +
  theme_bw()


## Use Gene Expression Data Now ##
load("~/Documents/Research/Druggable Genome/cmap/processed_data/X.RData")
drug.cols <- grep("_drug",colnames(X))

drug.seq <- c(1,10,25,50,100,200,500,1000)
n.reps <- 20
for(n.drug in drug.seq) {
  for(rep in n.reps) {
    beta <- rep(0,ncol(X))
    beta[-drug.cols] <- rnorm(n=(ncol(X)-length(drug.cols)),sd=0.5)
    pos.drugs <- sample(drug.cols,size=n.drug)
    beta[pos.drugs] <- rnorm(n=n.drug)
    y <- 1 + X %*% beta + rnorm(n=nrow(X),sd=0.25)
    fit <- glmnet(x=X,y=y,alpha=1,intercept=TRUE,standardize = TRUE)
    params <- lasic(fit,X,y,T)
  }
}




