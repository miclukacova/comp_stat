simulate <- function(n = 2000, par) {
  mu <- as.numeric(par[1])
  sigma2 <- as.numeric(par[2])
  nu <- as.numeric(par[3])
  w <- rchisq(n, df = nu)
  x <- rnorm(n, mean = mu, sd = sqrt(sigma2 * nu / w))
  return(list(w = w, x = x))
}


