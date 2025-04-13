#' Computes the autocorrelation function of a Latent autoregressive (LAR1) for quantiles of a LNAR1 process
#' @param l lag
#' @param a the (log) quantile of interest
#' @param phi correlation parameter of the latent AR(1)
#' @param sigma_e the standard deviation of the innovations of the latent AR(1)
#'
#' @return the autocorrelation at lag `l`
#' @export LAR_autocorr_2
#'
#' @details the LAR is a **binary** process.
#'  This is mostly an internal function.
#' Here we are computing the autocorrelation of a process of the type
#' exp(Z_t) <= log(a), where Z_t is an AR(1, lmeam, lsd, phi).
#' 
LAR_autocorr_2 <- function(l, a, phi, sigma_e){

  sigma_m <- sigma_e/sqrt(1 - phi^2)
  covk <-  phi^l * sigma_m^2
  Sigma <- matrix(c(sigma_m^2, covk, covk, sigma_m^2), nrow = 2)
  P1 <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                         upper = c(a, a),
                         mean = rep(0, 2),
                         sigma = Sigma)
  P2 <- mvtnorm::pmvnorm(lower = c(-Inf, a),
                         upper = c(a, Inf),
                         mean = rep(0, 2),
                         sigma = Sigma)
  ans <- (P1[1]/stats::pnorm(q = a, sd = sigma_m)) -
    (P2[1]/stats::pnorm(q = a, sd = sigma_m, lower.tail = FALSE))
  return(ans)
}
LAR_autocorr_2 <- Vectorize(LAR_autocorr_2)

#' Autocovariance function of LNAR1 process
#'
#' @param k the lag 
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#' @return the autocovariance at `k`
#' @export autocov_LNAR1
#'
autocov_LNAR1 <- function(k, lmean, lsd, phi){
  v <- lsd^2
  uy.sq <- exp(2*lmean + v)
  EW <- exp(2 * lmean + v * (1 + phi^k))
  ans <- EW - uy.sq
  return(ans)
}
autocov_LNAR1 <- Vectorize(autocov_LNAR1)

#' Autcorrelation function of a LNAR1 process
#'
#' @param k the lag
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#'
#' @return the autocorrelation at `k`
#' @export autocorr_LNAR1
#'
autocorr_LNAR1 <- function(k, lmean, lsd, phi){
  covar <- autocov_LNAR1(k, lmean, lsd, phi)
  v <- lsd^2
  vy <- (exp(v) - 1) * exp(2*lmean + v)
  ans <- covar/vy
  return(ans)
}
autocorr_LNAR1 <- Vectorize(autocorr_LNAR1)

#' Compute the efficiency (for the mean) of a LNAR1 process
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#' @param K an upper bound up to which to sum the autocorrelation function.
#' 
#' @export LNAR1_eff
#' @return the efficiency, which should be >= 0
#'
LNAR1_eff <- function(lmean, lsd, phi, K = 1e5){
  corrk <- function(k) autocorr_LNAR1(k = k,
                                       lmean = lmean, lsd = lsd, phi = phi)
  theoEffy <-  1 / (1 + 2 * sum(corrk(1:K)))
  return(theoEffy)
}

#' Find the autocorrelation parameter for a given efficiency
#'
#' @param eff the desired efficiency
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the value of phi (between -1 and 1)
#' @export find_phi_LNAR1 
#'
find_phi_LNAR1 <- function(eff, lmean, lsd){
  opt_eff <- function(cand){
    theo <- LNAR1_eff(lmean = lmean, lsd = lsd, phi = cand)
    return(
      (eff - theo)^2
    )
  }
  opt_eff <- Vectorize(opt_eff)
  Opt <- stats::optimise(opt_eff, interval = c(-1, 1))
  return(Opt$minimum)
}

#' Autocorrelation function for the p-th quantile of a LNAR1 process
#'
#' @param k lag
#' @param xi_p the *value* of the p-th quantile
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#' @param K an upper bound up to which to sum the autocorrelation function.
#'  Default is K = 5000
#' @return the autocorrelation value
#' @export LNAR_quantile_autocorr
#'
LNAR_quantile_autocorr <- function(k, xi_p, lmean, lsd, phi, K){
  
  v <- lsd^2 
  ## Depends on the LAR
  sgE <- sqrt(v) * sqrt(1 - phi ^ 2)
  
  g_j <- function(j){
    LAR_autocorr_2(l = j,
                   a = log(xi_p) - lmean,
                   sigma_e = sgE,
                   phi = phi)
  }
  g_j <- Vectorize(g_j)
  
  f_j <- function(j){
    return((1 - j/K) * g_j(j))
  }
  f_j <- Vectorize(f_j)
  
  return(f_j(k))
}

#' Compute the long-term (LT) variance of a quantile indicator
#'
#' @param p a quantile (between 0 and 1)
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#' @param K an upper bound up to which to sum the autocorrelation function.
#'  Default is K = 5000
#'
#' @return a scalar containt the LT variance
#' @export compute_LT_variance
#'
compute_LT_variance <- function(p, lmean, lsd, phi, K = 5000){
  x <- stats::qlnorm(p = p, 
                meanlog = lmean,
                sdlog = lsd)
  autcorrs <- LNAR_quantile_autocorr(k = 1:(K - 1),
                                     xi_p = x,
                                     lmean = lmean,
                                     lsd = lsd,
                                     phi = phi,
                                     K = K)
  
  SigmaSq_p <- p * (1 - p) * (1 + 2 * sum(autcorrs))
  
  return(SigmaSq_p)
}
#' Fits a log-normal autorregressive (LNAR) model of order 1
#'
#' @param X draws from MCMC
#'
#' @returns a list containing the estimated LNAR1 parameters (lmean, lsd and phi)
#' @export fit_LNAR1
#'
fit_LNAR1 <- function(X){
  fit <- stats::arima(log(X), order = c(1, 0, 0))
  
  m.hat <- fit$coef[2]
  phi.hat <- fit$coef[1]
  s.hat <- sqrt(fit$sigma2/(1-phi.hat^2))
  
  out <- list(
    mu_hat = as.numeric(m.hat),
    sd_hat = as.numeric(s.hat),
    phi_hat = as.numeric(phi.hat) 
  )
  return(out)
}