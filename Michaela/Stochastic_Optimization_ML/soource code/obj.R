##### Objective function #####################################

H <- function(par, x, y) {
  mean((y - f(par, x))^2)
} 

H_mult <- function(alpha, beta, gamma, rho, x, y){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) H(par, x, y))
}