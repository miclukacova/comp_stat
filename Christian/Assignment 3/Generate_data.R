parameters <- function(mu, sigma, nu){
  structure(
    list(
      par = c(mu, sigma, nu),
      mu = mu,
      sigma = sigma,
      nu = nu
    ),
    class = "parameters"
  )
}


simulate.parameters <- function(parameters_object, N = 1000, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  w <- rchisq(N, parameters_object$nu)
  x <- rnorm(N, mean = parameters_object$mu, sd = sqrt(parameters_object$nu * parameters_object$sigma^2 / w))
  
  return(list(w = w, x = x))
}


