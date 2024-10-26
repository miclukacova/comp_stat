mse_dens <- function(est, true) {
  mean((est - true)^2)
}
