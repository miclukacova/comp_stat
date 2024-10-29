
Qgrad <- function(par, nu, alpha_prime, beta_prime, data){
  # n <- length(x)
  # dQ_dmu <- - alpha_prime / (nu * sigma^2) * sum((x - mu) / beta_prime)
  # dQ_dsigma <- - n/sigma + alpha_prime / (nu * sigma^3) * sum( (x - mu)^2 / beta_prime)
  
  mu <- par[1]
  sigma <- par[2]
  return(-d_marginal_log_likelihood_cpp(x = data, mu = mu, sigma = sigma, nu = nu))
}



fisher_information_naive <- function(par, nu, alpha_prime, data){
  n <- length(x)
  mu <- par[1]
  sigma <- par[2]
  
  beta_prime <- beta_prime_func(x = x, nu = nu, mu = mu, sigma = sigma)
  grad <- -Qgrad(x = x, par = par, nu = nu, alpha_prime = alpha_prime, beta_prime = beta_prime)/n
  
  fisher_information <- matrix(0, nrow = 2, ncol = 2)
  
  for (i in seq_along(x)){
    
    fisher_i <- Qgrad(par, 
                      nu = nu, 
                      alpha_prime = alpha_prime, 
                      beta_prime = beta_prime[i],
                      x = data[i])
    
    fisher_information <- fisher_information + (fisher_i - grad) %*% t(fisher_i - grad)
  }
  
  return(fisher_information)
}



fisher_information_naive2 <- function(x, mu, sigma, nu){
  n <- length(x)
  n_1grad <- -d_marginal_log_likelihood_cpp(x = x, mu = mu, sigma = sigma, nu = nu)/n
  
  fisher_information <- matrix(0, nrow = 2, ncol = 2)
  for (i in seq_along(x)){
    
    fisher_i <- -d_marginal_log_likelihood_cpp(x = x[i], mu = mu, sigma = sigma, nu = nu)
    
    fisher_information <- fisher_information + (fisher_i - n_1grad) %*% t(fisher_i - n_1grad)
  }
  return(fisher_information)
}



fisher_information_naive3 <- function(x, mu, sigma, nu){
  fisher_information <- matrix(0, nrow = 2, ncol = 2)
  for (i in seq_along(x)){
    
    fisher_i <- -d_marginal_log_likelihood_cpp(x = x[i], mu = mu, sigma = sigma, nu = nu)
    
    fisher_information <- fisher_information + (fisher_i) %*% t(fisher_i)
  }
  return(fisher_information)
}


fisher_information <- function(x, mu, sigma, nu){
  
  single_fisher <- function(x, mu, sigma, nu){
    fisher_i <- -d_marginal_log_likelihood_cpp(x = x, mu = mu, sigma = sigma, nu = nu)
    fisher_i %*% t(fisher_i)
  }
  
  fisher_information_list <- lapply(x, function(x_i) single_fisher(x = x_i, mu = mu, sigma = sigma, nu = nu))
  Reduce('+', fisher_information_list)
}
