#### S3 class

parameters <- function(mu, sigma2, nu) {
  structure(
    list(
      mu = mu,
      sigma2 = sigma2,
      nu = nu, 
      pars = c(mu, sigma2, nu)),
    class = "My_parameters"
  )
}

sim <- function(x) {
  UseMethod("sim")
}

sim <- function(n, obj) {
  mu <- obj$mu
  sigma2 <- obj$sigma2
  nu <- obj$nu
  w <- rchisq(n, df = nu)
  x <- rnorm(n, mean = mu, sd = sqrt(sigma2 * nu / w))
  return(list(w = w, x = x))
}

mle <- function(x) {
  UseMethod("sim")
}

mle <- function(obj, data) {
  c(mle.mu(data$x, data$w), mle.sigma2(data$x, data$w, mle.mu(data$x, data$w), obj$nu))
}

print.My_parameters <- function(x) {
  cat("Parameters: \n")
  cat("mu: ", x$mu, "\n")
  cat("sigma2: ", x$sigma2, "\n")
  cat("nu: ", x$nu, "\n")
}
