# Loglikelihood

loglik <- function(x, par, nu) {
  mu <- par[1]
  sigma2 <- par[2]
  sum(- log(sqrt(sigma2)) - (nu + 1) / 2 * log(1 + (x - mu)^2 / (nu * sigma2)))
}

grad_loglik <- function(mu, sigma, nu, x){
  n <- length(x)
  d_mu <- (nu + 1) * sum((x - mu) / (nu * sigma^2 + (x - mu)^2))
  d_sigma <- - n/sigma + (nu + 1) / sigma * sum((x - mu)^2 / (nu * sigma^2 + (x - mu)^2))
  return(c(d_mu, d_sigma))
}

