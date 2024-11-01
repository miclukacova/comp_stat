## Fischer information

fis_inf_single_term <- function(x_i, mu, sigma, nu, full_x){
  N <- length(full_x)
  (grad_loglik(mu, sigma, nu, x_i) - 1 / N * grad_loglik(mu, sigma, nu, full_x)) %*%
    t(grad_loglik(mu, sigma, nu, x_i) - 1 / N * grad_loglik(mu, sigma, nu, full_x))
}

fis_inf <- function(mu, sigma, nu, full_x) {
  
  listt <- lapply(full_x, function(x_i) fis_inf_single_term(x_i, mu, sigma, nu, full_x))
  
  Reduce("+", listt)
}

grad_Q <- function(mu, sigma, nu, x){
  n <- length(x)
  E.W <- E.step(x, c(mu, sigma^2), nu)
  d_mu <- sum(E.W * (x - mu) / (nu * sigma^2))
  d_sigma <- - n/sigma + sum(E.W * (x - mu)^2 / (nu * sigma^3))
  return(c(d_mu, d_sigma))
}

