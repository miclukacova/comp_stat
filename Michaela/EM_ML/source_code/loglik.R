# Loglikelihood

loglik <- function(x, par, nu) {
  mu <- par[1]
  sigma2 <- par[2]
  sum(- log(sqrt(sigma2)) - (nu + 1) / 2 * log(1 + (x - mu)^2 / (nu * sigma2)))
}

