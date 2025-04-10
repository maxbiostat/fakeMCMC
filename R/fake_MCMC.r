#' Generate 'fake' MCMC draws from a latent gaussian autorregressive order 1 (LNAR1) process
#'
#' @param N number of draws
#' @param lmean log mean
#' @param lsd log standard deviation
#' @param phi correlation parameter of the latent AR(1)
#' @param LN whether to generate from a log-normal (default is true)
#'
#' @returns a vector of size N with the draws
#' @export generate_fake_MCMC
#'
generate_fake_MCMC <- function(N = 1E4,
                               lmean = 0,
                               lsd = 1,
                               phi = -0.75,
                               LN = TRUE) {
  v <- lsd^2
  sigma <- sqrt(v)
  sigma_e <- lsd * sqrt(1 - phi ^ 2)
  rawSimu <- stats::arima.sim(n = N, list(ar = c(phi)), sd = sigma_e) + lmean
  if (LN) {
    exampleSimu <- exp(rawSimu)
  } else{
    exampleSimu <- rawSimu
  }
  return(exampleSimu)
}