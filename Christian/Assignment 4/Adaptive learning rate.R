decay_scheduler <- function(gamma0 = 1, a = 1, K = 1, gamma1, n1)
{
  force(a)
  if (!missing(gamma1) && !missing(n1))
    K <- n1^a * gamma1 / (gamma0 - gamma1)
  b <- gamma0 * K
  function(n) b / (K + n^a)
}



batch <- function(
    par,
    samp,
    gamma,
    grad,          # Function of parameter and observation
    m = 50,        # Mini-batch size
    ...
){
  M <- floor(length(samp) / m) 
  for (j in 0:(M - 1)) {
    i <- samp[(j * m + 1):(j * m + m)]
    par <- par - gamma * grad(par, i, ...)
  }
  par
}



momentum <- function() {
  rho <- 0
  function(
    par,
    samp,
    gamma,
    grad,
    m = 50,
    beta = 0.95,
    ...
    # Mini-batch size
    # Momentum memory
  ){
    M <- floor(length(samp) / m) 
    for (j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      # Using '<<-' assigns the value to rho in the enclosing
      environment
      rho <<- beta * rho + (1 - beta) * grad(par, i, ...)
      par <- par - gamma * rho
    }
    par
  } 
}



adam <- function() {
  rho <- v <- 0
  function(
    par,
    samp,
    gamma,
    grad,
    m = 50,
    beta1 = 0.9,
    beta2 = 0.9,
    ...
    # Mini-batch size
    # Momentum memory
    # Second moment memory
  ){
    M <- floor(length(samp) / m) 
    for (j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      gr <- grad(par, i, ...)
      rho <<- beta1 * rho + (1 - beta1) * gr
      v <<- beta2 * v + (1 - beta2) * gr^2
      par <- par - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par
  } 
}