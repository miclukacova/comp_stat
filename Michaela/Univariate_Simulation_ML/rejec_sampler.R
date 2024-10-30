# Rejection sampler

gauss_rej <- function(mu, sigma, alpha, f) {
  sampler <-function(n){
    y <- numeric(n)
    counter <- 0
    for (i in 1:n) {
      reject <- TRUE
      while (reject) {
        y0 <- rnorm(1, mu, sigma)
        u <- runif(1)
        reject <- u > f(y0) / dnorm(y0, mu, sigma) * alpha
        counter <- counter + 1
      }
      y[i] <- y0
    }
    return(list("sample" = y, "accept" = n / counter))
  }
  return(sampler)
}


#gauss_rej <- function(n, mu, sigma, alpha, f) {
#  y <- numeric(n)
#  counter <- 0
#  for (i in 1:n) {
#    reject <- TRUE
#    while (reject) {
#      y0 <- rnorm(1, mu, sigma)
#      u <- runif(1)
#      reject <- u > f(y0) / dnorm(y0, mu, sigma) * alpha
#      counter <- counter + 1
#    }
#    y[i] <- y0
#  }
#  return(list("sample" = y, "accept" = n / counter))
#}#