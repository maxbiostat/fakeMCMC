#' Compute the mean of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the mean
#' @export lognormal_mean
#'
#' @examples 
#' lognormal_mean(0, 1)
#' mean(exp(rnorm(1e6)))
lognormal_mean <- function(lmean, lsd){
  exp(lmean + lsd^2/2) 
}
#' Compute the median of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the median
#' @export lognormal_median
#'
#' @examples
#' lognormal_median(0, 1)
#' median(exp(rnorm(1e6)))
lognormal_median <- function(lmean, lsd){
  exp(lmean) 
}
#' Compute the variance of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the variance
#' @export lognormal_variance
#'
#' @examples 
#' lognormal_variance(0, 1)
#' var(exp(rnorm(1e6)))
lognormal_variance <- function(lmean, lsd){
  (exp(lsd^2) - 1) * exp(2*lmean + lsd^2) 
}
#' Compute the HPD of a log-normal distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return a vector of size two with the HPD endpoints
#' @export lognormal_hpd
#'
#' @examples
#' lognormal_hpd(0.95, 0, 1)
#' lognormal_hpd(0.90, 0, 1)
lognormal_hpd <- function(alpha = 0.95, lmean, lsd) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- stats::qlnorm(p = c(q1, q2),
                   meanlog = lmean,
                   sdlog = lsd)
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-20,
                  upper = 1 - alpha,
                  tol = 1E-30)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(qlnorm(
    p = c(q1.opt, q2.opt),
    meanlog = lmean,
    sdlog = lsd
  ))
}

#' Get the percentiles corresponding to the HPD of a log-normal
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return
#' @export
#'
#' @examples
#' lognormal_hpd_percentiles(0.95, 0, 1)
#' lognormal_hpd_percentiles(0.90, 0, 1)
#' lognormal_hpd_percentiles(0.95, 0, .1)
#' lognormal_hpd_percentiles(0.90, 0, .1)
lognormal_hpd_percentiles <- function(alpha = 0.95, lmean, lsd){
  hpd <- lognormal_hpd(alpha = alpha, lmean, lsd)
  return(
    stats::plnorm(q = hpd, meanlog = lmean, sdlog = lsd)
  )
}

#' Compute the skewness of a log-normal
#'
#' @param lsd log standard deviation
#'
#' @return a scalar with the skewness
#' @export lognormal_skewness
#'
#' @examples
#' lognormal_skewness(.1)
#' lognormal_skewness(.5)
#' lognormal_skewness(1)
#' curve(lognormal_skewness)
lognormal_skewness <- function(lsd){
  v <- lsd^2
  skw <- (exp(v) + 2) * sqrt(exp(v) - 1)
  return(skw)
}
lognormal_skewness <- Vectorize(lognormal_skewness)
#' Mode of a log-normal
#'
#' @param lmean log-normal mean parameter (mu)
#' @param lsd log-normal standard deviation parameter (sigma)
#'
#' @return the mode
#' @export lognormal_mode
#'
lognormal_mode <- function(lmean = 0, lsd = 1) {
  
  if (lsd <= 0)
    stop("lsd must be > 0")
  
  mo <- exp(lmean - lsd^2)
  
  return(mo)
}
#' Compute the central BCI for a log-normal distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return a vector of size two with the BCI endpoints
#' @export lognormal_bci
#'
#' @examples
#' lognormal_bci(0.95, 0, 1)
#' lognormal_bci(0.90, 0, 1)
lognormal_bci <- function(alpha, lmean, lsd) {
  ps <- c((1 - alpha), (1 + alpha))/2
  return(stats::qlnorm(p = ps, meanlog = lmean, sdlog = lsd))
}

#'  Compute the HPD of an exponential distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param theta the rate parameter.
#'
#' @return a vector of size two with the HPD endpoints
#' @export exponential_hpd
#'
#' @examples
#' exponential_hpd(0.95, 1)
#' exponential_hpd(0.90, 1)
#' exponential_hpd(0.95, 10)
#' exponential_hpd(0.90, 10)
exponential_hpd <- function(alpha = 0.95, theta) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- stats::qexp(p = c(q1, q2), rate = theta)
    return(cand[2] - cand[1])
  }
  Opt <- stats::optimise(opt_int,
                  a = alpha,
                  lower = 1E-20,
                  upper = 1 - alpha,
                  tol = 1E-30)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(qexp(
    p = c(q1.opt, q2.opt), rate = theta))
}
#'  Compute the central BCI of an exponential distribution
#'
#' @param alpha the BCI level  (default is 0.95)
#' @param theta the rate parameter.
#'
#' @return a vector of size two with the HPD endpoints
#' @export exponential_bci
#'
#' @examples
#' exponential_bci(0.95, 1)
#' exponential_bci(0.90, 1)
#' exponential_bci(0.95, 10)
#' exponential_bci(0.90, 10)
exponential_bci <- function(alpha = 0.95, theta) {
  if(theta <= 0) stop("Theta must be strictly positive")
  ps <- c((1 - alpha), (1 + alpha))/2
  return(stats::qexp(p = ps, rate = theta))
}