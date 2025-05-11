#' Estimate the parameters of a LNAR process
#'
#' @param X draws from the observed process
#' @param level confidence level for parameter confidence intervals
#'
#' @return a data.frame with point estimates and confidence intervals for all
#'  three parameters
#' @export est_LNAR_pars
#'
est_LNAR_pars <- function(X, level = 0.95) {
  
  M <- length(X)
  Y <- log(X)
  fit <- stats::arima(Y, order = c(1, 0, 0))
  coefs <- c(stats::coef(fit), stats::sd(Y))
  ints <- stats::confint(fit)
  x2qs <- stats::qchisq(p = c(1 - level, 1 + level)/2, df = M - 1)
  S2 <- stats::var(Y)
  out <- data.frame(
    parameter = c("phi", "mu", "sigma"),
    point = coefs,
    lwr = c(ints[, 1], sqrt((M - 1)/x2qs[2] * S2)),
    upr = c(ints[, 2], sqrt((M - 1)/x2qs[1] * S2))
  )
  rownames(out) <- NULL
  return(out)
}