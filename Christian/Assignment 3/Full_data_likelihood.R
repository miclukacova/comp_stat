

full_data_log_likelihood <- function(x, w, mu, sigma, nu){
  n <- length(x)
  
  term1 <- -n/2 * log(pi * nu * sigma^2)
  term2 <- - n * (nu + 1)/2 * log(2)
  term3 <- - n * lgamma(nu/2)
  term4 <- sum( (nu - 1)/2 * log(w) - w / 2 * (1 + (x - mu)^2 / (nu * sigma^2)) )
  
  log_likelihood <- term1 + term2 + term3 + term4
  return(log_likelihood)
}


full_data_mle <- function(x, w, nu){
  
  n <- length(x)
  
  mu_hat <- sum(w * x) / sum(w)
  sigma_hat <- sqrt(sum(w * (x - mu_hat)^2) / (n * nu))
  
  return(list(mu = mu_hat, sigma = sigma_hat))
  
}


FDLL_function_factory <- function(x, w, nu){
  
  n <- length(x)
  
  # Negative full data log-likelihood
  H <- function(par){
    - full_data_log_likelihood(x = x, w = w, mu = par[[1]], sigma = par[[2]], nu = nu) / n
  }
  
  # Negative gradient of marginal log-likelihood
  grad <- function(par){
    - d_marginal_log_likelihood(x = x, par[[1]], par[[2]], nu = nu) / n
  }
  
  return(list(H = H, grad = grad))
}
