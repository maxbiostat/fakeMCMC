#' Computes the autocorrelation function of a Latent autoregressive (LAR1) for quantiles of a LNAR1 process
#' @param l lag
#' @param a the (log) quantile of interest
#' @param phi correlation parameter of the latent AR(1)
#' @param sigma_e the standard deviation of the innovations of the latent AR(1)
#'
#' @return the autocorrelation at lag `l`
#' @export LAR_autocorr
#'
#' @details the LAR is a **binary** process.
#'  This is mostly an internal function.
LAR_autocorr <- function(l, a, phi, sigma_e){
  
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
LAR_autocorr <- Vectorize(LAR_autocorr)

#' Find the 'alpha' parameter of LAR(1) process
#'
#' @param sigma_e  the standard deviation of the innovations of the latent AR(1)
#' @param phi correlation parameter of the latent AR(1)
#' @param p the marginal probability that X = 1
#'
#' @return a scalar between 0 and 1
#' @export compute_LAR_alpha
#' @details
#' This is the approximation of a LAR(1) by a first-order two-state Markov chain
#' 
#'
compute_LAR_alpha <- function(sigma_e, phi, p){
  sigma_marginal <- sigma_e/sqrt(1 - phi^2)
  a <- stats::qnorm(p = p, sd = sigma_marginal)
  ff <- function(z){
    exp(stats::dnorm(z, mean = 0, sd = sigma_marginal, log = TRUE) + 
          stats::pnorm(q = a, mean = phi*z, sd = sigma_e, log.p = TRUE) )
  }
  num <- stats::integrate(ff, a, Inf)$value
  ldenom <- stats::pnorm(a, sd = sigma_marginal,
                         lower.tail = FALSE, log.p = TRUE)
  lans <- log(num) - ldenom
  return(exp(lans))
}

#' Find the 'beta' parameter of a LAR(1) process
#'
#' @param sigma_e the standard deviation of the innovations of the latent AR(1)
#' @param phi correlation parameter of the latent AR(1)
#' @param p the marginal probability that X = 1
#'
#' @return a scalar between 0 and 1
#'
#' @export compute_LAR_beta
#' @details
#'  This is the approximation of a LAR(1) by a first-order two-state Markov chain
#' 
compute_LAR_beta <- function(sigma_e, phi, p){
  sigma_marginal <- sigma_e/sqrt(1-phi^2)
  a <- stats::qnorm(p = p, sd = sigma_marginal)
  gg <- function(z){
    exp(stats::dnorm(z, mean = 0, sd = sigma_marginal, log = TRUE) + 
          stats::pnorm(q = a, mean = phi*z, sd = sigma_e,
                 log.p = TRUE, lower.tail = FALSE) )
  }
  num <- stats::integrate(gg, -Inf, a)$value
  ldenom <- stats::pnorm(a, sd = sigma_marginal, log.p = TRUE)
  lans <- log(num) - ldenom
  return(exp(lans))
}

#' Compute the true effective sample size (ESS) for an LAR(1) process
#'
#' @param N the number of draws
#' @param p the marginal probability that X = 1
#' @param phi correlation parameter of the latent AR(1)
#' @param sigma_e the standard deviation of the innovations of the latent AR(1)
#' @param K  an upper bound up to which to sum the autocorrelation function.
#'  Default is K = 5000
#'
#' @return a scalar with the ESS
#' @export LAR_ess
#'
LAR_ess <-  function(N, p, phi, sigma_e, K = 1000){
  sigma_marginal <- sigma_e/sqrt(1-phi^2)
  S <- sum(LAR_autocorr(l = 1:K,
                        a = stats::qnorm(p = p, sd = sigma_marginal),
                        phi = phi, sigma_e = sigma_e))
  return(N/(1 + 2*S))
}

#' Sample from a LAR(1) process
#'
#' @param N the number of draws
#' @param phi correlation parameter of the latent AR(1)
#' @param target.p the marginal probability that X = 1
#'
#' @return a vector with \code{N} draws
#' @export generate_LAR_TS
#'
#' @examples
#' p <- .23
#' thePhi <- get_LAR_phi(Eff = .6, p = p)
#' X <- generate_LAR_TS(N = 1E3, phi = thePhi, target.p = p)
generate_LAR_TS <- function(N, phi, target.p){
  Z <- stats::arima.sim(n = N,
                        list(ar = c(phi)),
                 sd = sqrt((1 - phi^2))) # latent process
  X <- as.vector( (Z <= stats::qnorm(p = target.p)) + 0 )
  return(X)
}

#' Computes the phi associated with a given pair {ESS, p}
#'
#' @param Eff the target efficiency
#' @param p the marginal probability that X = 1
#'
#' @return a number between -1 and 1 
#' @export get_LAR_phi
#' @details
#' Assumes sigma_e^2 = 1-rho^2, meaning the marginal is N(0, 1)
#' 
get_LAR_phi <- function(Eff, p){
  
  obj <- function(rho){
    tent_eff <- LAR_ess(N = 1, p = p,
                        phi = rho, sigma_e = sqrt(1 - rho^2))
    return((tent_eff - Eff)^2)
  }
  obj <- Vectorize(obj)
  Opt <- stats::optimise(obj, interval = c(-.999, .999))
  return(Opt$minimum)
}

#' Generate draws from a latent autorregressive (LAR1) model
#'
#' @param N number of draws
#' @param phi correlation parameter between -1 and 1
#' @param target.p target marginal probability that X = 1
#'
#' @return a vector of \code{N} binary draws
#' @export generate_LAR_TS
#'
#' @examples
#' the.phi <- get_LAR_phi(Eff = .2, p = .12)
#' X <- generate_LAR_TS(N = 1e3, phi = the.phi, target.p = .12)
#' acf(X)
generate_LAR_TS <- function(N, phi, target.p){
  Z <- stats::arima.sim(n = N, list(ar = c(phi)),
                 sd = sqrt((1 - phi^2))) # latent process
  X <- as.vector( (Z <= stats::qnorm(p = target.p)) + 0 )
  return(X)
}