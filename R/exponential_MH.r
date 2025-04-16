#' Independence Metropolis-Hastings with exponential(1) target
#'
#' @param M number of steps
#' @param theta proposal rate parameter
#' @param x0 initial state
#'
#' @returns a data.frame containing the states, the log target density and whether the move was accepted
#' @export MH_exponential
#'
MH_exponential <- function(M = 1e4,
                           theta = 1,
                           x0 = 0.1){
  ## Jones & Hobert (2001) Stat Sci
  
  res <- data.frame(matrix(NA, ncol = 3, nrow = M + 1))
  
  curr <- x0
  curr_ltarget <- stats::dexp(x = curr, rate = 1, log = TRUE)
  names(res) <- c("x", "lp", "accp")
  
  res[1, ] <- c(curr, curr_ltarget, NA)
  
  for(i in 2:(M+1)){
    
    prop <- stats::rexp(1, rate = theta)
    prop_ltarget <- stats::dexp(x = prop, rate = 1, log = TRUE)
    
    log_accp <- (prop_ltarget + stats::dexp(curr, rate = theta, log = TRUE)) -
                (curr_ltarget + stats::dexp(prop, rate = theta, log = TRUE))
  
    if(log(stats::runif(1)) < log_accp){
      ## Accept
      res[i, ] <- c(prop, prop_ltarget, 1)
      curr_ltarget <- prop_ltarget
      curr <- prop
    }else{
      ## Reject
      res[i, ] <- c(curr, curr_ltarget, 0)
    }
  }
  return(res)
}