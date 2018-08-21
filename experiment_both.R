require(MASS)
mu.true <- c(0, 1)
# mu1=0, mu2=1
Sigma.true <- matrix(c(1, 1.6, 1.6, 4), nrow = 2)
# sigma11=1, sigma22=4, rho=0.8

theta <- function(mu, Sigma) {
  c(mu, Sigma[-3], Sigma[2] / sqrt(Sigma[1] * Sigma[4]))
}


theta.true <- theta(mu.true, Sigma.true)

n <- 100
B <- 1000

Bucks_Method <- function(sample.missing.both) {
  complete_cases <- complete.cases(sample.missing.both)
  y1_miss <- is.na(sample.missing.both[, 1])
  y2_miss <- is.na(sample.missing.both[, 2])
  mean.complete <- colMeans(sample.missing.both[complete_cases, ])
  Sigma.complete <- cov(sample.missing.both[complete_cases, ])
  beta_21 <- Sigma.complete[2] / Sigma.complete[1]
  beta_12 <- Sigma.complete[2] / Sigma.complete[4]
  sample.imputed <- sample.missing.both
  # \hat y1 = \bar y1_com +beta1 (y2 - \bar y2_com)
  sample.imputed[y1_miss, 1] <- mean.complete[1] +
    beta_12 * (sample.missing.both[y1_miss, 2] - mean.complete[2])
  sample.imputed[y2_miss, 2] <- mean.complete[2] +
    beta_21 * (sample.missing.both[y2_miss, 1] - mean.complete[1])
  theta(colMeans(sample.imputed), cov(sample.imputed))
}

EM_algorithm <- function(sample.missing.both, eps = 1e-7) {
  n <- nrow(sample.missing.both)
  complete_cases <- complete.cases(sample.missing.both)
  y1_miss <- is.na(sample.missing.both[, 1])
  y2_miss <- is.na(sample.missing.both[, 2])
  mean.new <- colMeans(sample.missing.both[complete_cases, ])
  cov.new <- cov(sample.missing.both[complete_cases, ])[c(1, 2, 4)]
  mean.old <- rep(Inf, 2)
  cov.old <- rep(Inf, 3)
  #s的固定部分
  s_1.com <- sum(sample.missing.both[complete_cases, 1])
  s_2.com <- sum(sample.missing.both[complete_cases, 2])
  s_11.com <- sum(sample.missing.both[complete_cases, 1] ^ 2)
  s_22.com <- sum(sample.missing.both[complete_cases, 2] ^ 2)
  s_12.com <- sum(sample.missing.both[complete_cases, 1] *
                    sample.missing.both[complete_cases, 2])
  while (max(abs(c(mean.new - mean.old, cov.new - cov.old))) > eps) {
    mean.old <- mean.new
    cov.old <- cov.new
    #先计算y1观测而y2缺失的案例
    beta_21_1 <- cov.old[2] / cov.old[1]
    beta_20_1 <- mean.old[2] - beta_21_1 * mean.old[1]
    sigma_22_1 <- cov.old[3] - cov.old[2] ^ 2 / cov.old[1]
    s_1 <- s_1.com + sum(sample.missing.both[y2_miss, 1])
    s_11 <- s_11.com + sum(sample.missing.both[y2_miss, 1] ^ 2)
    hat_y2 <-
      beta_20_1 + beta_21_1 * sample.missing.both[y2_miss, 1]
    s_2 <- s_2.com + sum(hat_y2)
    s_22 <- s_22.com + sum(hat_y2 ^ 2 + sigma_22_1)
    s_12 <- s_12.com + sum(hat_y2 * sample.missing.both[y2_miss, 1])
    #再计算y2观测而y1缺失的案例
    beta_11_2 <- cov.old[2] / cov.old[3]
    beta_10_2 <- mean.old[1] - beta_11_2 * mean.old[2]
    sigma_11_2 <- cov.old[1] - cov.old[2] ^ 2 / cov.old[3]
    s_2 <- s_2 + sum(sample.missing.both[y1_miss, 2])
    s_22 <- s_22 + sum(sample.missing.both[y1_miss, 2] ^ 2)
    hat_y1 <-
      beta_10_2 + beta_11_2 * sample.missing.both[y1_miss, 2]
    s_1 <- s_1 + sum(hat_y1)
    s_11 <- s_11 + sum(hat_y1 ^ 2 + sigma_11_2)
    s_12 <- s_12 + sum(hat_y1 * sample.missing.both[y1_miss, 2])
    # 重新估计参数
    mean.new <- c(s_1, s_2) / n
    cov.new <- c(s_11 / n - mean.new[1] ^ 2,
                 s_12 / n - mean.new[1] * mean.new[2],
                 s_22 / n - mean.new[2] ^ 2)
  }
  c(mean.new, cov.new, cov.new[2] / sqrt(cov.new[1] * cov.new[3]))
}





# Bucks_Method(sample.missing.both)
# theta(colMeans(sample.missing.both[complete_cases,]),
#       cov(sample.missing.both[complete_cases,]))
# EM_algorithm(sample.missing.both)

both <-
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
    
    sample.missing.both <- sample.complete
    if (mechanism == "MAR") {
      #前50个中y2较大的30个
      y1.NA.cases <-
        head(order(sample.complete[1:(n / 2), 2], decreasing = T), n / 10 * 3)
      #后50个中y1较小的30个
      y2.NA.cases <-
        head(order(sample.complete[(n / 2 + 1):n, 1], decreasing = T), n / 10 *
               3)+n/2
    }
    else{
      # MCAR
      y1.NA.cases <- 1:(3 / 10 * n) # 前30个
      y2.NA.cases <- (n - 3 / 10 * n + 1):n # 后30个
    }
    sample.missing.both[y1.NA.cases, 1] <- NA
    sample.missing.both[y2.NA.cases, 2] <- NA
    complete_cases <- complete.cases(sample.missing.both)
    
    #均值填补
    sample.mean_imputed <- sample.missing.both
    sample.mean_imputed[y1.NA.cases, 1] <-
      mean(sample.complete[-y1.NA.cases, 1])
    sample.mean_imputed[y2.NA.cases, 2] <-
      mean(sample.complete[-y2.NA.cases, 2])
    theta.mean_imputed <- theta(colMeans(sample.mean_imputed),
                                cov(sample.mean_imputed))
    
    rbind(
      theta.complete = theta.complete,
      theta.CC = theta(colMeans(sample.missing.both[complete_cases,]),
                       cov(sample.missing.both[complete_cases,])),
      theta.mean_imputed = theta.mean_imputed,
      theta.MLE = EM_algorithm(sample.missing.both),
      theta.Bucks = Bucks_Method(sample.missing.both)
    )
  }

#both(mechanism = "MCAR")
ptm<-proc.time()
estimate<-array(dim=c(B,5,6))
for (i in 1:B){
  estimate[i,,]= both(mechanism = "MAR",
                          seed=i)
}
proc.time()-ptm

estimate.mean<-apply (estimate,c(2,3),mean)
colnames(estimate.mean)<-c("mu1","mu2","sigma_11","sigma_12",
                           "sigma_22","rho")
rownames(estimate.mean)<-c("complete","CC",
                           "uncond mean","EM","Buck's")
estimate.sd<-apply (estimate,c(2,3),sd)
colnames(estimate.sd)<-c("mu1","mu2","sigma_11","sigma_12",
                         "sigma_22","rho")
rownames(estimate.sd)<-c("complete","CC",
                         "uncond mean","EM","Buck's")

estimate.mean
estimate.sd