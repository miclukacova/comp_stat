grad <- function(par, x, y) {
  
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  const_term <- exp(beta * log(x) - alpha)
  
  d_alpha <- (rho - gamma) / (1 + const_term)^2 * const_term
  d_beta <- - log(x) * (rho - gamma) / (1 + const_term)^2 * const_term
  d_gamma <- 1 - 1 / (1 + const_term)
  d_rho <- 1 / (1 + const_term)
  
  return(- 2 * c(d_alpha, d_beta, d_gamma, d_rho) * (y - f(par, x)))
}

grad_mult <- function(par, x, y) {
  rowSums(sapply(seq_along(x), function(i) grad(par, x[i], y[i])))
}
