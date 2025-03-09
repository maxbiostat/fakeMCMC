generate_fake_MCMC <- function(N = 1E4,
                               m = 0,
                               v = 1,
                               phi = -0.75,
                               LN = TRUE) {
  sigma <- sqrt(v)
  sigma_e <- sigma * sqrt(1 - phi ^ 2)
  rawSimu <- stats::arima.sim(n = N, list(ar = c(phi)), sd = sigma_e) + m
  if (LN) {
    exampleSimu <- exp(rawSimu)
  } else{
    exampleSimu <- rawSimu
  }
  return(exampleSimu)
}