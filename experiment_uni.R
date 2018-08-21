require(MASS)
mu.true <- c(0, 1)
# mu1=0, mu2=1
Sigma.true <- matrix(c(1, 1.6, 1.6, 4), nrow = 2)
# sigma11=1, sigma22=4, rho=0.8

theta <- function(mu, Sigma) {
  c(mu, Sigma[-3], Sigma[2] / sqrt(Sigma[1] * Sigma[4]))
}


theta.true <- theta(mu.true, Sigma.true)

n<-100
B<-1000

MLE.monotne <- function(sample.missing.uni) {
  complete_cases <- complete.cases(sample.missing.uni)
  mu1 <- mean(sample.missing.uni[, 1])
  sigma_11 <- var(sample.missing.uni[, 1])
  mean_obs <- colMeans(sample.missing.uni[complete_cases, ])
  Sigma_obs <- cov(sample.missing.uni[complete_cases, ])
  beta <- Sigma_obs[2] / Sigma_obs[1]
  sigma_12 <- beta * sigma_11
  mu2 <- mean_obs[2] + beta * (mu1 - mean_obs[1])
  sigma_22 <- Sigma_obs[4] + beta ^ 2 * (sigma_11 - Sigma_obs[1])
  c(mu1,
    mu2,
    sigma_11,
    sigma_12,
    sigma_22,
    sigma_12 / sqrt(sigma_11 * sigma_22))
}

regression.impute<-function(sample.missing.uni) {
  complete_cases <- complete.cases(sample.missing.uni)
  sample<-data.frame(y1=sample.missing.uni[,1],
                     y2=sample.missing.uni[,2])
  lm.obs<-lm(y2~y1,data=sample,subset=complete_cases)
  sample[-which(complete_cases),2]<-
    predict(lm.obs,newdata=data.frame(
      y1=sample[-which(complete_cases),1]))
  theta(colMeans(sample),cov(sample))
}

monotone <-
  function(n = 100,
           mu = mu.true,
           Sigma = Sigma.true,
           mechanism = "MAR",
           seed = 1) {
    set.seed(seed)
    sample.complete <-
      mvrnorm(n = n, mu = mu, Sigma = Sigma)
    theta.complete <- theta(colMeans(sample.complete),
                            cov(sample.complete))
    
    set.seed(seed)
    if(mechanism=="MAR") {
      NA.cases <- head(order(sample.complete[, 1], decreasing = T), n / 2)
    }
    else{
      NA.cases <- sample(n, n / 2)
    }
                       
    sample.missing.uni <- sample.complete
    sample.missing.uni[NA.cases, 2] <- NA
    
    
    sample.mean_imputed <- sample.missing.uni
    sample.mean_imputed[NA.cases, 2] <-
      mean(sample.complete[-NA.cases, 2])
    theta.mean_imputed <- theta(colMeans(sample.mean_imputed),
                                cov(sample.mean_imputed))
    
    rbind(
      theta.complete = theta.complete,
      theta.CC = theta(colMeans(sample.missing.uni[-NA.cases, ]),
                       cov(sample.missing.uni[-NA.cases, ])),
      theta.mean_imputed = theta.mean_imputed,
      theta.MLE = MLE.monotne(sample.missing.uni),
      reg_impute= regression.impute(sample.missing.uni)
    )
  }

estimate<-array(dim=c(B,5,6))
for (i in 1:B){
  estimate[i,,]= monotone(mechanism = "MAR",
                          seed=i)
}

estimate.mean<-apply (estimate,c(2,3),mean)
colnames(estimate.mean)<-c("mu1","mu2","sigma_11","sigma_12",
                           "sigma_22","rho")
rownames(estimate.mean)<-c("complete","CC",
                           "uncond mean","MLE","reg imputed")
estimate.sd<-apply (estimate,c(2,3),sd)
colnames(estimate.sd)<-c("mu1","mu2","sigma_11","sigma_12",
                           "sigma_22","rho")
rownames(estimate.sd)<-c("complete","CC",
                           "uncond mean","MLE","reg imputed")

estimate.mean
estimate.sd


