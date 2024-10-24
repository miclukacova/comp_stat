## Sampler  ##############################################

gauss_sample <- function(N, par, omega = 1) {
  log_x <- rnorm(N, 0, omega^2)
  x <- exp(log_x)
  y <- f(par = par, x = x) + rnorm(N, 0, 0.5)
  return(data.frame(x = x, y = y))
}

grid_sample <- function(N, par) {
  grid <- exp(1:15)
  x <- sample(grid, N, replace = TRUE)
  y <- f(par = par, x = x) + rnorm(N, 0, 0.5)
  return(data.frame(x = x, y = y))
}

###### Parameters class ###########################################

parameters <- function(alpha, beta, gamma, rho) {
  structure(
    list(
      alpha = alpha,
      beta = beta,
      gamma = gamma, 
      rho = rho,
      par = c(alpha, beta, gamma, rho)),
    class = "My_params"
  )
}

sim <- function(x) {
  UseMethod("sim")
}

sim <- function(object, N, omega = 1, grid = FALSE, scale = FALSE) {
  if(grid){
    data <- grid_sample(N = N, par = object$par)
  }
  data <- gauss_sample(N = N, par = object$par, omega)
  #For scaling
  if(scale){
    scaled_data <- scale(data) 
    # Shift the data so that the minimum value is 1 (or any positive value)
    data <- scaled_data + abs(min(scaled_data)) + 1
  }
  return(data)
}
