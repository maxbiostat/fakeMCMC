#' Run a fake Metropolis (Fabreti-Hoehna) algorithm with fixed acceptance probability and a given
#' autocorrelation time
#'
#' @param M  number of steps
#' @param rgen the function which generates random variables from the target
#' @param theta the parameters which control rgen, e.g. theta = c(0, 1) for a standard normal
#' @param tau the autocorrelation time, a positive scalar.
#'
#' @returns a vector with M samples from the fake Metropolis sampler.
#' @export fake_Metrop
#' @examples
#' gaussian_gen <- function(theta) rnorm(1, mean = theta[1], sd = theta[2])
#' Xs <- fake_Metrop(M = 100, rgen = gaussian_gen, theta = c(0, 1), tau = 3)
#' hist(Xs)
fake_Metrop <- function(M, rgen, theta, tau){
  # alpha <- 2*target.ess/(N + target.ess)
  alpha <- 2/(1 + tau)
  log_alpha <- log(alpha)
  X <- rep(NA, M)
  X[1] <-  rgen(theta)
  for (i in 2:M){
    u <- runif(1)
    if(log_alpha >= log(stats::runif(1)) ){
      X[i] <- rgen(theta)
    }else{
      X[i] <- X[i-1]
    }
  }
  return(X)
}