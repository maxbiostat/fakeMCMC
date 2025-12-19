#' Exponential autorregressive model EAR(1)
#' Lewis, P. A. (1978). Simple models for positive-valued and discrete-valued time series with  ARMA correlation structure.
#' 
#' @param M number of steps
#' @param beta acceptance probability( a scalar between 0 and 1)
#' @param x0 initial state
#'
#' @returns a data.frame containing the states, the log target density and whether the move was accepted
#' @export run_EAR1
#'
run_EAR1 <- function(M = 1e4, beta, x0 = 0.1){
 
  if(beta < 0 || beta >= 1) stop("beta needs to be in [0, 1)")
  
  log_beta <- log(beta)
  
  res <- data.frame(matrix(NA, ncol = 3, nrow = M + 1))
  names(res) <- c("x", "lp", "accp")
  
  res[1, ] <- c(x0, NA, NA)
  
  for(i in 2:(M+1)){
    
    if(log(stats::runif(1)) < log_beta){
      
      res[i, ] <- c(beta * res[i-1, 1], NA, 1)
      
    }else{
      
      res[i, ] <- c(beta * res[i-1, 1] + stats::rexp(1), NA, 0)
    }
    
  }
  res[ , 2] <- stats::dexp(x = res[, 1], rate = 1, log = TRUE)
  return(res)
}