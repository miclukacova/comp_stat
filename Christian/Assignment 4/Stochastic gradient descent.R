log_logistic_dose_response_model <- function(x, par){
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  return(gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha)))
}



gradient <- function(par, i, x, y,...){
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  x_i <- x[i]
  y_i <- y[i]
  
  expbetalogxalpha <- exp(beta * log(x_i) - alpha)
  
  identical_part <- - 2 * (y_i - log_logistic_dose_response_model(x_i, par))
  
  grad_alpha <- mean(identical_part * (rho - gamma) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_beta <- - mean(identical_part * (rho - gamma) * log(x[i]) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_gamma <- mean(identical_part * (1 - 1 / (1 + expbetalogxalpha)))
  grad_rho <- mean(identical_part / (1 + expbetalogxalpha))
  
  return(c(grad_alpha, grad_beta, grad_gamma, grad_rho))
}




sgd <- function(
    par,
    N, # Sample size
    gamma, # Decay schedule or a fixed learning rate
    epoch = NULL,
    maxit = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    ...) {
  
  if (is.function(gamma)) gamma <- gamma(1:maxit) 
  gamma <- rep_len(gamma, maxit)
  
  for (n in 1:maxit) {
    if (!is.null(cb)) cb()
    
    samp <- sampler(N)
    
    if (is.null(epoch)){
      for (j in 1:N) {
        i <- samp[j]
        par <- par - gamma[n] * grad(par, i, ...)
      }
    } else {
      par <- epoch(par, samp, gamma[n], ...)
    }
  }
  par
}



