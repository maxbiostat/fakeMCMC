my_mcse.q <- 
function (x, q, size = NULL, g = NULL,
          method = c("bm", "obm", "sub"), warn = FALSE) 
{
  method = match.arg(method)
  if (is.function(g)) 
    x = sapply(x, g)
  n = length(x)
  if (n < 1000) {
    if (warn) 
      warning("too few samples (less than 1,000)")
    if (n < 10) 
      return(NA)
  }
  if (is.null(size)) {
    b <- batchSize(x = x, method = method)
  }
  else if (size == "sqroot") {
    b <- floor(sqrt(n))
  }
  else if (size == "cuberoot") {
    b <- floor(n^(1/3))
  }
  else {
    if (!is.numeric(size) || size < 1 || size >= n || floor(n/size) <= 
        1) {
      warning("size is either too large, too small, or not a number. Setting 'size' to n^(1/2)")
      size = sqrt(n)
    }
    b <- floor(size)
  }
  a <- floor(n/b)
  if (!is.numeric(q) || q <= 0 || q >= 1) {
    stop("'q' must be from (0, 1)")
  }
  xi.hat = quant(x, q)
  if (method == "bm") {
    var.hat = mcseqbm(x, b, xi.hat)
    lx <- log(x)
    lnpars <- c(mean(lx), sd(lx))
    f.hat = dlnorm(x = xi.hat, meanlog = lnpars[1], sdlog = lnpars[2])
    se = sqrt(var.hat/n)/f.hat
  }
  else if (method == "obm") {
    var.hat = mcseqobm(x, b, xi.hat)
    lx <- log(x)
    lnpars <- c(mean(lx), sd(lx))
    f.hat = dlnorm(x = xi.hat, meanlog = lnpars[1], sdlog = lnpars[2])
    se = sqrt(var.hat/n)/f.hat
  }
  else {
    var.hat = mcseqsub(x, b, q, quant)
    se = sqrt(var.hat/n)
  }
  value = list(est = xi.hat, se = se, nsim = n)
  value
}
environment(my_mcse.q) <- asNamespace('mcmcse')