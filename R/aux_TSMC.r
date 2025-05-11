#' Generate draws from a two state Markov chain
#'
#' @param N the number of draws
#' @param target.ess target effective sample size (ESS)
#' @param target.p the marginal probability that X = 1
#' @param x0 initial state. Either 0, 1 or \code{NULL}.
#' if \code{NULL}, a fair coin is flipped to pick \code{x0}
#'
#' @returns a vector with \code{N} 
#' @export generate_MC_TS
#'
#' @examples
#' X <- generate_MC_TS(N = 1e3, target.ess = 150, target.p = .23)
generate_MC_TS <- function(N, target.ess, target.p = 0.5, x0 = NULL) {
  requireNamespace("BinaryMarkovChains")
  target.act <- N / target.ess
  target.alpha <- (2 * target.p) / (target.act + 1)
  ## The line below clamps out attempts to set efficiencies higher than
  ### theory allows
  maxA <- BinaryMarkovChains::get_max_alpha(target.p)
  target.alpha <- min(target.alpha, maxA)
  target.beta <-
    exp(log(target.alpha) + log(1 - target.p) - log(target.p))
  target.mat <- matrix(
    c(1 - target.alpha, target.alpha,
      target.beta, 1 - target.beta),
    ncol = 2,
    nrow = 2,
    byrow = TRUE
  )
  
  target.TS.MC <- new(
    "markovchain",
    states = c("0", "1"),
    transitionMatrix = target.mat,
    name = "Binary"
  )
  if(is.null(x0)) x0 <- stats::rbinom(1, 1, .5)
  X <- as.numeric(markovchain::markovchainSequence(
    n = N,
    markovchain = target.TS.MC,
    t0 = paste(x0)
  ))
  
  return(X)
}