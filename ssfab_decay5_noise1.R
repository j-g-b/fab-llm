#
source("helpers.R")

#
n1 <- 50
n2 <- 50
d <- 10
p <- 200
Rs <- c(0.8, 0.66, 0.5)

nsim <- 100
nsigs <- 20

#
mu1 <- rep(0, p)
E <- svd(matrix(rnorm(p*d), nrow = p))[["u"]]

#
Deltas <- list()
Pwrs <- list()
for(k in 1:length(Rs)){
  Delta <- matrix(NA, p, nsigs)
  Pwr <- matrix(NA, 3, nsigs)
  r <- Rs[k]
  for(i in 1:nsigs){
    #
    Sigma_d <- diag(seq(1, 0.1, length.out = p)^5)
    Sigma_d <- 50*(Sigma_d / sum(Sigma_d))
    Sigma_U <- svd(matrix(rnorm(p*p), nrow = p))[["u"]]
    mu2 <- E%*%rnorm(d) + rnorm(p)
    mu2 <- mu2 / sqrt(sum(mu2^2))
    if(i > 15){
      mu2 <- rep(0, p)
    }
    Delta[, i] <- mu2
    
    p_rp <- rep(NA, nsim)
    p_fab <- rep(NA, nsim)
    p_ssfab <- rep(NA, nsim)
    for(j in 1:nsim){
      #
      X1 <- mu1%*%t(rep(1, n1)) + Sigma_U%*%sqrt(Sigma_d)%*%matrix(rnorm(p*n1), nrow = p)
      X2 <- mu2%*%t(rep(1, n2)) + Sigma_U%*%sqrt(Sigma_d)%*%matrix(rnorm(p*n2), nrow = p)
      
      #
      Z <- matrix(rnorm(p*d), nrow = p)
      ZX1 <- crossprod(Z, X1)
      ZX2 <- crossprod(Z, X2)
      zxbar1 <- rowMeans(ZX1)
      zxbar2 <- rowMeans(ZX2)
      SZ <- ((n1-1)*cov(t(ZX1)) + (n2-1)*cov(t(ZX2)))/(n1 + n2 - 2)
      
      #
      RP_stat <- ((n1+n2-d-1)/(d*(n1+n2-2)))*(n1*n2/(n1+n2))*crossprod(zxbar2 - zxbar1, crossprod(solve(SZ), zxbar2 - zxbar1))
      p_rp[j] <- 1 - pf(RP_stat, d, n1 + n2 - 1 - d)
      
      #
      EX1 <- crossprod(E, X1)
      EX2 <- crossprod(E, X2)
      Exbar1 <- rowMeans(EX1)
      Exbar2 <- rowMeans(EX2)
      SE <- ((n1-1)*cov(t(EX1)) + (n2-1)*cov(t(EX2)))/(n1 + n2 - 2)
      
      #
      FAB_stat <- ((n1+n2-d-1)/(d*(n1+n2-2)))*(n1*n2/(n1+n2))*crossprod(Exbar2 - Exbar1, crossprod(solve(SE), Exbar2 - Exbar1))
      p_fab[j] <- 1 - pf(FAB_stat, d, n1 + n2 - 1 - d)
      
      #
      ssn1 <- sample(n1, floor(r*n1))
      ssn2 <- sample(n2, floor(r*n2))
      
      nssX1 <- X1[, setdiff(1:n1, ssn1)]
      nssX2 <- X2[, setdiff(1:n2, ssn2)]
      Stilde <- ((n1 - floor(r*n1)-1)*cov(t(nssX1)) + (n2 - floor(r*n2)-1)*cov(t(nssX2)))/(n1 - floor(r*n1) + n2 - floor(r*n1) - 2)
      Stildei <- MASS::ginv(Stilde + (1e-1)*(p/(n1 - floor(r*n1) + n2 - floor(r*n2)))*diag(1, p))
      
      ssEX1 <- crossprod(Stildei%*%E, X1[, ssn1])
      ssEX2 <- crossprod(Stildei%*%E, X2[, ssn2])
      ssExbar1 <- rowMeans(ssEX1)
      ssExbar2 <- rowMeans(ssEX2)
      ssSE <- ((floor(r*n1)-1)*cov(t(ssEX1)) + (floor(r*n2)-1)*cov(t(ssEX2)))/(floor(r*n1) + floor(r*n1) - 2)
      
      #
      ssFAB_stat <- ((floor(r*n1)+floor(r*n2)-d-1)/(d*(floor(r*n1)+floor(r*n2)-2)))*(floor(r*n1)*floor(r*n2)/(floor(r*n1)+floor(r*n2)))*crossprod(ssExbar2 - ssExbar1, crossprod(solve(ssSE), ssExbar2 - ssExbar1))
      p_ssfab[j] <- 1 - pf(ssFAB_stat, d, floor(r*n1) + floor(r*n2) - 1 - d)
      
    }
    Pwr[1, i] <- mean(p_rp < 0.05)
    Pwr[2, i] <- mean(p_fab < 0.05)
    Pwr[3, i] <- mean(p_ssfab < 0.05)
    print(i)
  }
  Deltas[[k]] <- Delta
  Pwrs[[k]] <- Pwr
}

#
saveRDS(list(Deltas, Pwrs), "results/ssfab_decay5_noise1.rds")
