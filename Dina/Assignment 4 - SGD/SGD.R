log_logistic_dose_response_model <- function(x, alpha, beta, gamma, rho){
  return(gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha)))
}


gradient <- function(par, x, y){
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  expbetalogxalpha <- exp(beta * log(x) - alpha)
  
  identical_part <- - 2 * (y - log_logistic_dose_response_model(x, alpha, beta, gamma, rho))
  
  grad_alpha <- identical_part * (rho - gamma) * expbetalogxalpha / (1 + expbetalogxalpha)^2
  grad_beta <- - identical_part * (rho - gamma) * log(x) * expbetalogxalpha / (1 + expbetalogxalpha)^2
  grad_gamma <- identical_part * (1 - 1 / (1 + expbetalogxalpha))
  grad_rho <- identical_part / (1 + expbetalogxalpha)
  
  return(c(grad_alpha, grad_beta, grad_gamma, grad_rho))
}

sgd <- function(
    par,
    grad, # Function of parameter and observation index
    data, # Data
    n, # Sample size
    gamma, # Decay schedule or a fixed learning rate
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    ...) {
  x <- data$x
  y <- data$y
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter) {
    if (!is.null(cb)) cb()
    samp <- sampler(n)
    for (j in 1:n) {
      i <- samp[j]
      x_i <- x[i]
      y_i <- y[i]
      par <- par - gamma[k] * grad(par, x_i, y_i)
    }
  }
  par
}