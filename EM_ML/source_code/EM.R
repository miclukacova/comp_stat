## EM algorithm ##################################################

## E step

E.W <- function(X, par){
  (nu + 1) / (1 + (X - par[1])^2 / (nu * par[2]))
}

#E.logW <- function(x, parameters){
#  mu <- parameters$mu
#  sigma2 <- parameters$sigma2
#  nu <- parameters$nu
#  digamma ((nu + 1) / 2) + log(2) - log(1 + (x - mu)^2 / (nu * sigma2))
#}
#
#Q <- function(X, parameters){
#  mu <- parameters$mu
#  sigma2 <- parameters$sigma2
#  nu <- parameters$nu
#  
#  first_term <- - n * length(X) * log(sqrt(pi * nu * sigma2) * 2^((nu + 1) / 2) * gamma(nu / 2))
#  second_term <- sum(E.logW(X, parameters)) * (nu - 1) / 2
#  third_term <- -sum(E.W(X, parameters) / 2 * (1 + (X - mu)^2 / (nu * sigma2)))
#  
#  first_term + second_term + third_term
#}

## M step

M.step <- function(E.W, X) {
  mu_k <- sum(E.W * X) / sum(E.W)
  sigma2_k <- sum(E.W * (X - mu_k)^2) / (length(X) * nu)
  return(c(mu_k, sigma2_k))
}


## EM algorithm

em_factory <- function(e_step, m_step, eps = 1e-6) {
  force(e_step); force(m_step); force(eps)
  function(par, X, epsilon = eps, cb = NULL, ...) {
    repeat {
      par0 <- par
      par <- m_step(E.W = e_step(X = X, par = par), X = X)
      if (!is.null(cb)) cb()
      if (sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon))
        break
    }
    par  # Returns the parameter estimate
  }
}
