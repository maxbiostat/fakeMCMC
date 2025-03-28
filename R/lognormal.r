#' Compute the BCI for a log-normal
#'
#' @param alpha a number between 0 and 1 corresponding to the level of the BCI
#' @param lmean log-normal mean parameter (mu)
#' @param lsd log-normal standard deviation parameter (sigma)
#'
#' @return a vector of size 2 with the BCI
#' @export lognormal_bci
#'
lognormal_bci <- function(alpha = 0.95,
                          lmean = 0,
                          lsd = 1) {
  if (alpha < 0 || alpha > 1)
    stop("alpha must be in (0, 1) ")
  
  if (lsd <= 0)
    stop("lsd must be > 0")
  
  return(
    stats::qlnorm(
      p = c(1 - alpha, 1 + alpha) / 2,
      meanlog = lmean,
      sdlog = lsd
    )
  )
}

#' Compute the highest density interval for a log-normal
#'
#' @param alpha a number between 0 and 1 corresponding to the level of the BCI
#' @param lmean log-normal mean parameter (mu)
#' @param lsd log-normal standard deviation parameter (sigma)
#'
#' @return a vector of size 2 with the HDI/HPD
#' @export lognormal_hpd
#'
lognormal_hpd <- function(alpha = 0.95,
                          lmean = 0,
                          lsd = 1) {
  if (alpha < 0 || alpha > 1)
    stop("alpha must be in (0, 1) ")
  
  if (lsd <= 0)
    stop("lsd must be > 0")
  
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- stats::qlnorm(p = c(q1, q2),
                   meanlog = lmean,
                   sdlog = lsd)
    return(cand[2] - cand[1])
  }
  Opt <- stats::optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(
    stats::qlnorm(
      p = c(q1.opt, q2.opt),
      meanlog = lmean,
      sdlog = lsd
    )
  )
}

#' Compute the quantiles corresponding the the HPD endpoints in a log-normal
#'
#'
#' @param alpha a number between 0 and 1 corresponding to the level of the BCI
#' @param lmean log-normal mean parameter (mu)
#' @param lsd log-normal standard deviation parameter (sigma)
#'
#' @return a vector of size 2 with the quantiles
#' @export lognormal_hpd_percentiles
#'
lognormal_hpd_percentiles <-
  function(alpha = 0.95,
           lmean = 0,
           lsd = 1) {
    if (alpha < 0 || alpha > 1)
      stop("alpha must be in (0, 1) ")
    
    if (lsd <= 0)
      stop("lsd must be > 0")
    
    hpd <- lognormal_hpd(alpha = alpha, lmean, lsd)
    return(
      stats::plnorm(
        q = hpd,
        meanlog = lmean,
        sdlog = lsd
      )
    )
  }

#' Mean (expectation) of a log-normal
#'
#' @param lmean log-normal mean parameter (mu)
#' @param lsd log-normal standard deviation parameter (sigma)
#'
#' @return the mean
#' @export lognormal_mean
#'
lognormal_mean <- function(lmean = 0, lsd = 1) {
  
  if (lsd <= 0)
    stop("lsd must be > 0")
  
  m <- exp(lmean + lsd^2/2)
  return(m)
}

#' Skewness of a log-normal distribution
#'
#' @param lsd  log-normal standard deviation parameter (sigma)
#'
#' @return a scalar with the skewness
#' @export lognormal_skewness
#'
lognormal_skewness <- function(lsd = 1) {
  if (lsd <= 0)
    stop("lsd must be > 0")
  
  v <- lsd ^ 2
  skw <- (exp(v) + 2) * sqrt(exp(v) - 1)
  return(skw)
}
lognormal_skewness <- Vectorize(lognormal_skewness)

#' Median of a log-normal distribution
#'
#' @param lmean  log-normal mean parameter (mu)
#'
#' @return a scalar with the median
#' @export lognormal_median
#'
lognormal_median <- function(lmean = 0) {
  md <- exp(lmean)
  return(md)
}
lognormal_median <- Vectorize(lognormal_median)

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
