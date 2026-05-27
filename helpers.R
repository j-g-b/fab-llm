#
msqrt <- function(M, inv = F){
  Msvd <- svd(M)
  if(inv){
    Msvd[["u"]]%*%diag(1/sqrt(Msvd[["d"]]))%*%t(Msvd[["u"]])
  } else {
    Msvd[["u"]]%*%diag(sqrt(Msvd[["d"]]))%*%t(Msvd[["u"]])
  }
}

#
T_fab <- function(y, X, Psi, ss0, null = F, S = 1000){
  #
  n <- length(y)
  PsiXinv <- solve(solve(Psi) + crossprod(X)/ss0)
  
  #
  if(!null){
    u <- y / sqrt(sum(y^2))
    x <- c(sqrt(1/ss0 - crossprod(u, X)%*%PsiXinv%*%crossprod(X, u)/ss0^2))
    return(-n*log(x))
  } else {
    Y <- matrix(rnorm(n*S), nrow = S)
    U <- t(apply(Y, 1, function(y){y / sqrt(sum(y^2))}))
    xU <- c(apply(U, 1, function(u){sqrt(1/ss0 - crossprod(u, X)%*%PsiXinv%*%crossprod(X, u)/ss0^2)}))
    return(-n*log(xU))
  }
  
}

#
T_F <- function(y, X){
  n <- nrow(X)
  p <- ncol(X)
  IPy <- y - X%*%(MASS::ginv(crossprod(X))%*%crossprod(X, y))
  c(((n-p)/p)*(t(y)%*%y - t(y)%*%IPy)/(t(y)%*%IPy))
}

#
#
min_lam <- function(z, p){
  optim(0.5, function(lambda){(z - p*log(lambda))/(1 - lambda)}, method = "Brent", lower = 0, upper = 1)$par
}

#
pr_bound <- function(z, p){
  ml <- min_lam(z, p)
  pgamma((z - p*log(ml))/(1 - ml), p/2, 1/2, lower.tail = F)
}

#
c_alpha <- function(alpha, p){
  zu <- 10*p
  zl <- -10*p
  prb <- pr_bound((zu + zl)/2, p)
  while(abs(prb - alpha) > 1e-4){
    if(prb > alpha){
      zl <- (zu + zl)/2
      prb <- pr_bound((zu + zl)/2, p)
    } else {
      zu <- (zu + zl)/2
      prb <- pr_bound((zu + zl)/2, p)
    }
  }
  return((zu + zl)/2)
}

#
find_signal_noise_mle <- function(Y, U, V, R, sigmasq){
  require(dfoptim)
  # Get dimensions
  n <- nrow(U)
  m <- nrow(V)
  p <- ncol(V)
  q <- ncol(U)
  # Mean-center Y
  Y <- Y - mean(Y)
  # Compute SVD of U, V
  normY <- sum(Y^2)
  if(!all(dim(U) == c(1,1))){
    svdU <- svd(U)
  }
  if(!all(dim(V) == c(1,1))){
    svdV <- svd(V)
  }
  #
  if(!all(dim(U) == c(1,1))){
    QY <- t(svdU[["u"]])%*%Y
  } else {
    QY <- Y
  }
  #
  if(!all(dim(V) == c(1,1))){
    QY <- QY%*%svdV[["u"]]
  }
  QY <- c(QY)
  normQY <- sum(QY^2)
  #
  if(!all(dim(V) == c(1,1))){
    svs <- svdV[["d"]]
    if(m > p){
      svs <- c(svs, rep(0, m - p))
    }
  } else {
    svs <- 1
  }
  #
  if(!all(dim(U) == c(1,1))){
    usvs <- svdU[["d"]]
    if(n > q){
      usvs <- c(usvs, rep(0, n - q))
    }
    svs <- c(sapply(svs, function(x){x*usvs}))
  }
  #
  Rmean <- 1/mean(1/R)
  N <- prod(dim(Y))
  non_zero_svs <- svs != 0
  any_zero_svs <- any(!non_zero_svs)
  sigmasq <- sigmasq / Inf
  #
  opt <- function(taupsi){
    tausq <- taupsi[1]
    psisq <- taupsi[2]
    (0.5)*sum(log(((psisq)*(svs^2) + (sigmasq / Rmean) + tausq))) + # log-determinant
      (0.5)*sum((1 / ((psisq)*(svs[non_zero_svs]^2) + (sigmasq / Rmean) + tausq))*(QY^2)) + # quad form directions of X
      (0.5)*ifelse(any_zero_svs, 1, 0)*(1 / ((sigmasq / Rmean) + tausq))*(normY - normQY) +# quad form anti directions of X
      (0.95)*(normY / N)/tausq +
      (0.05)*(normY / sum(svs^2))/psisq
  }
  #
  opt_res <- dfoptim::nmkb(c(0.1, 0.001), opt, lower = c(0, 0))
  opt_tausq <- opt_res[["par"]][1]
  opt_psisq <- opt_res[["par"]][2]
  return(list(tausq = opt_tausq, psisq = opt_psisq))
}


Rcpp::sourceCpp("fabp_lin_reg.cpp")
#
#
#' Compute FAB p-values under linear regression linking model
#'
#' @description Comutes FAB p-values based on a linear regression linking model for matrix data given row and column covariates.
#'
#' @param Y matrix of averages of experimental readout values over R replicates
#' @param S matrix of standard errors of experimental readout values over R replicates
#' @param R number of replicates per readout value
#' @param U the row features
#' @param V the column features
#' @param pool_sampling_var logical; indicates whether sampling variance should be assumed the same across hypothesis tests
#' @param Y1 matrix of averages of experimental readout values over R1 replicates (used for contrast scores for two-sample t-tests; Y is taken to be the other sample)
#' @param S1 matrix of standard errors of experimental readout values in Y1 over R1 replicates
#' @param R1 number of replicates per readout value in Y1
#'
#' @return A data.frame of FAB p-values and the standard UMP p-values, one for each entry in Y (or each contrast score Y - Y1)
#'
#' @export fabp_lin_reg
#'
fabp_lin_reg <- function(Y, S, R, U, V, pool_sampling_var = T, Y1 = NULL, S1 = NULL, R1 = NULL){
  # Get experimental dimensions
  vecY <- c(Y)
  vecR <- c(R)
  vecS <- c(S)
  m <- nrow(V)
  n <- nrow(U)
  p <- ncol(V)
  q <- ncol(U)
  # Split data and estimate sigmasq_hat and sigmasq_tild
  hat_indices <- sample(m*n, round(m*n / 2))
  tild_indices <- setdiff(1:(m*n), hat_indices)
  if(is.null(Y1)){
    hat_nu <- sum(vecR[hat_indices] - 1)
    tild_nu <- sum(vecR[tild_indices] - 1)
    sigmasq_hat <- sum((vecR[hat_indices] - 1)*vecS[hat_indices]^2) / hat_nu
    sigmasq_tild <- sum((vecR[tild_indices] - 1)*vecS[tild_indices]^2) / tild_nu
  } else {
    hat_nu <- sum(vecR[hat_indices] - 1) + sum(c(R1)[hat_indices] - 1)
    tild_nu <- sum(vecR[tild_indices] - 1) + sum(c(R1)[tild_indices] - 1)
    sigmasq_hat <- sum((vecR[hat_indices] - 1)*(vecS[hat_indices]^2) + (c(R1)[hat_indices] - 1)*((c(S1)[hat_indices])^2)) / hat_nu
    sigmasq_tild <- sum((vecR[tild_indices] - 1)*(vecS[tild_indices]^2) + (c(R1)[tild_indices] - 1)*((c(S1)[tild_indices])^2)) / tild_nu
    vecY <- vecY - c(Y1)
    Y <- Y - Y1
    R <- 1/((1/R) + (1/R1))
  }
  # Split data to estimate signal to noise
  larger_dim <- ifelse(nrow(Y) >= ncol(Y), "row", "col")
  if(larger_dim == "row"){
    #
    rand_rows <- sample(1:nrow(Y), round(nrow(Y)/2))
    mle_taupsi <- find_signal_noise_mle(Y[rand_rows, , drop = F], U[rand_rows, ], V, R[rand_rows, , drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_rows <- setdiff(1:nrow(Y), rand_rows)
    mle_taupsi <- find_signal_noise_mle(Y[not_rand_rows, , drop = F], U[not_rand_rows, ], V, R[not_rand_rows, , drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    not_rand_index <- c(sapply(not_rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  } else {
    #
    rand_cols <- sample(1:ncol(Y), round(ncol(Y)/2))
    mle_taupsi <- find_signal_noise_mle(Y[, rand_cols, drop = F], U, V[rand_cols, ], R[, rand_cols, drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_cols <- setdiff(1:ncol(Y), rand_cols)
    mle_taupsi <- find_signal_noise_mle(Y[, not_rand_cols, drop = F], U, V[not_rand_cols, ], R[, not_rand_cols, drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    not_rand_index <- c(sapply(not_rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  }
  #
  opt_tausq <- c(opt_tausq1, opt_tausq2)
  opt_psisq <- c(opt_psisq1, opt_psisq2)
  #
  if(!is.null(Y1)){
    vR <- 1/((1/vecR) + (1/c(R1)))
  } else {
    vR <- vecR
  }
  linking_estimators <- rcpp_fabp_lin_reg(vecY, sigmasq_tild, opt_tausq, opt_psisq, vR, U, V, PermIndx = perm_indx)
  # Extract linking model estimators
  theta_hat <- linking_estimators[[2]]
  tau_hat <- linking_estimators[[1]]
  # Compute LOOCV and LOOR^2
  loocv <- mean((vecY - theta_hat)^2)
  loor2 <- cor(vecY, theta_hat)^2
  # Compute test t-statistics and corresponding UMP, FAB p-values
  if(is.null(Y1)){
    if(pool_sampling_var){
      tstat <- vecY/(sqrt(sigmasq_hat/vR))
    } else {
      tstat <- vecY/(sqrt(vecS^2/vR))
    }
  } else {
    if(pool_sampling_var){
      tstat <- vecY/sqrt(sigmasq_hat/vR)
    } else {
      vecSsq <- ((vecR-1)*(vecS^2) + c(R1-1)*(c(S1)^2))/(vecR + c(R1) - 2)
      tstat <- vecY/sqrt(vecSsq/vR)
    }
  }
  b <- 2*theta_hat*sqrt(sigmasq_tild/vR)/tau_hat
  if(is.null(Y1)){
    if(pool_sampling_var){
      p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu)))
      fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu)))
    } else {
      p_values <- (1 - abs(pt(tstat, df = vecR - 1) - pt(-tstat, df = vecR - 1)))
      fabp_values <- (1 - abs(pt(tstat + b, df = vecR - 1) - pt(-tstat, df = vecR - 1)))
    }
  } else {
    if(pool_sampling_var){
      p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu)))
      fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu)))
    } else {
      p_values <- (1 - abs(pt(tstat, df = vecR + c(R1) - 2) - pt(-tstat, df = vecR + c(R1) - 2)))
      fabp_values <- (1 - abs(pt(tstat + b, df = vecR + c(R1) - 2) - pt(-tstat, df = vecR + c(R1) - 2)))
    }
  }
  #
  return(data.frame(row = rep(row.names(Y), m), column = rep(colnames(Y), each = n),
                    observed = vecY, predicted = theta_hat, model_error = mean(linking_estimators[[3]]), theta_var = tau_hat, error_var = mean(c(sigmasq_hat, sigmasq_tild)),
                    statistic = tstat, guess = b,
                    p = p_values, fabp = fabp_values,
                    fdr_p = p.adjust(p_values, method = "BH"), fdr_fabp = p.adjust(fabp_values, method = "BH"),
                    mse = loocv, r2 = loor2))
}

#
#
#' Compute FAB p-values under linear regression linking model
#'
#' @description Comutes FAB p-values based on a linear regression linking model for matrix data given row and column covariates.
#'
#' @param Y matrix of averages of experimental readout values over R replicates
#' @param S matrix of standard errors of experimental readout values over R replicates
#' @param R number of replicates per readout value
#' @param U the row features
#' @param V the column features
#' @param pool_sampling_var logical; indicates whether sampling variance should be assumed the same across hypothesis tests
#' @return A data.frame of FAB p-values and the standard UMP p-values, one for each entry in Y (or each contrast score Y - Y1)
#'
#' @export fabp_lin_reg
#'
fabp_lin_reg_z <- function(Y, S, R, U, V, pool_sampling_var = T){
  # Get experimental dimensions
  vecY <- c(Y)
  vecR <- c(R)
  vecS <- c(S)
  m <- nrow(V)
  n <- nrow(U)
  p <- ncol(V)
  q <- ncol(U)
  # Split data and estimate sigmasq_hat and sigmasq_tild
  hat_indices <- sample(m*n, round(m*n / 2))
  tild_indices <- setdiff(1:(m*n), hat_indices)
  hat_nu <- length(hat_indices)
  tild_nu <- length(tild_indices)
  sigmasq_hat <- sum(vecS[hat_indices]^2) / hat_nu
  sigmasq_tild <- sum(vecS[tild_indices]^2) / tild_nu
  
  # Split data to estimate signal to noise
  larger_dim <- ifelse(nrow(Y) >= ncol(Y), "row", "col")
  if(larger_dim == "row"){
    #
    rand_rows <- sample(1:nrow(Y), round(nrow(Y)/2))
    mle_taupsi <- find_signal_noise_mle(Y[rand_rows, , drop = F], U[rand_rows, ], V, R[rand_rows, , drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_rows <- setdiff(1:nrow(Y), rand_rows)
    mle_taupsi <- find_signal_noise_mle(Y[not_rand_rows, , drop = F], U[not_rand_rows, ], V, R[not_rand_rows, , drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    not_rand_index <- c(sapply(not_rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  } else {
    #
    rand_cols <- sample(1:ncol(Y), round(ncol(Y)/2))
    mle_taupsi <- find_signal_noise_mle(Y[, rand_cols, drop = F], U, V[rand_cols, ], R[, rand_cols, drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_cols <- setdiff(1:ncol(Y), rand_cols)
    mle_taupsi <- find_signal_noise_mle(Y[, not_rand_cols, drop = F], U, V[not_rand_cols, ], R[, not_rand_cols, drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    not_rand_index <- c(sapply(not_rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  }
  #
  opt_tausq <- c(opt_tausq1, opt_tausq2)
  opt_psisq <- c(opt_psisq1, opt_psisq2)
  #
  
  vR <- vecR
  linking_estimators <- rcpp_fabp_lin_reg(vecY, sigmasq_tild, opt_tausq, opt_psisq, vR, U, V, PermIndx = perm_indx)
  
  # Extract linking model estimators
  theta_hat <- linking_estimators[[2]]
  tau_hat <- linking_estimators[[1]]
  # Compute LOOCV and LOOR^2
  loocv <- mean((vecY - theta_hat)^2)
  loor2 <- cor(vecY, theta_hat)^2
  # Compute test t-statistics and corresponding UMP, FAB p-values
  if(pool_sampling_var){
    tstat <- vecY/(sqrt(sigmasq_hat/vR))
  } else {
    tstat <- vecY/(sqrt(vecS^2/vR))
  }
  b <- 2*theta_hat*sqrt(sigmasq_tild/vR)/tau_hat
  if(pool_sampling_var){
    p_values <- (1 - abs(pnorm(tstat) - pnorm(-tstat)))
    fabp_values <- (1 - abs(pnorm(tstat + b) - pnorm(-tstat)))
  } else {
    p_values <- (1 - abs(pnorm(tstat) - pnorm(-tstat)))
    fabp_values <- (1 - abs(pnorm(tstat + b) - pnorm(-tstat)))
  }
  
  #
  return(data.frame(row = rep(row.names(Y), m), column = rep(colnames(Y), each = n),
                    observed = vecY, predicted = theta_hat, model_error = mean(linking_estimators[[3]]), theta_var = tau_hat, error_var = mean(c(sigmasq_hat, sigmasq_tild)),
                    statistic = tstat, guess = b,
                    p = p_values, fabp = fabp_values,
                    fdr_p = p.adjust(p_values, method = "BH"), fdr_fabp = p.adjust(fabp_values, method = "BH"),
                    mse = loocv, r2 = loor2))
}
