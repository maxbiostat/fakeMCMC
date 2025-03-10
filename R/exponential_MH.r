MH_exponential <- function(M = 1e4,
                           theta = 1,
                           x0 = 0.1){
  ## TODO find a good citation
  
  res <- data.frame(matrix(NA, ncol = 3, nrow = M))
  
  curr <- x0
  curr_ltarget <- dexp(x = curr, rate = 1, log = TRUE)
  names(res) <- c("x", "lp", "accp")
  
  res[1, ] <- c(curr, curr_ltarget, 1)
  
  for(i in 2:M){
    
    prop <- rexp(1, rate = theta)
    prop_ltarget <- dexp(x = prop, rate = 1, log = TRUE)
    
    log_accp <- (prop_ltarget + dexp(curr, rate = theta, log = TRUE)) -
                (curr_ltarget + dexp(prop, rate = theta, log = TRUE))
  
    if(log(runif(1)) < log_accp){
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