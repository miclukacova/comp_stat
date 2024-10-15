##################################################
################### Packages #####################
##################################################

library(bench)
library(ggplot2)
library(tidyverse)
library(testthat)
library(profvis)
library(gridExtra)
theme_set(theme_bw())
library(CSwR)

## SGD ###############################################

sgd <- function(
    par,
    grad, # Function of parameter and observation index
    n, # Sample size
    gamma, # Decay schedule or a fixed learning rate
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    ...) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter) {
    if (!is.null(cb)) cb()
    samp <- sampler(n)
    for (j in 1:n) {
      i <- samp[j]
      par <- par - gamma[k] * grad(par, i, ...)
    }
  }
  par
}

## Function  ###########################################

f <- function(par, x) {
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  gamma + (rho - gamma)/(1 + exp(beta * log(x) - alpha))
}

## Gradient  ###########################################

grad <- function(par, i) {
  
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  d_alpha <- (rho - gamma)/(1 + exp(beta * log(x[i]) - alpha))^2 * exp(beta * log(x[i]) - alpha)
  d_beta <- - log(x[i]) * (rho - gamma)/(1 + exp(beta * log(x[i]) - alpha))^2 *
    exp(beta * log(x[i]) - alpha)
  d_gamma <- 1 - 1/(1 + exp(beta * log(x[i]) - alpha))
  d_rho <- 1/(1 + exp(beta * log(x[i]) - alpha))

  return(- c(d_alpha, d_beta, d_gamma, d_rho) * (y[i] - f(par, x[i])))
}

## Sampler  ##############################################

gauss_sample <- function(N, par, omega = 3) {
  log_x <- rnorm(N, 0, omega^2)
  x <- exp(log_x)
  y <- f(par, x) + rnorm(N, 0, 1)
  return(data.frame(x = x, y = y))
}

grid_sample <- function(N, par) {
  grid <- seq(exp(1), exp(10), length = 10^4)
  sample(grid, N, replace = TRUE)
}

## Gradient Descent #####################################

grad_desc <- function(
    par,
    grad,
    H,
    t0 = 1e-2,
    maxit = 1000,
    cb = NULL,
    epsilon = 1e-4,
    beta = 0.8,
    alpha = 0.1,
    ...) {
  for (i in 1:maxit) {
    
    # Calculations of objective and gradient
    value <- H(par)
    gr <- grad(par)
    grad_norm <- sum(gr^2)
    
    # Callback
    if (!is.null(cb)) cb()
    
    # Convergence criterion based on gradient norm
    if (grad_norm <= epsilon) break
    t <- t0
    # Proposed descent step
    par_new <- par - t * gr
    
    # Backtracking line search
    while (H(par_new) > value - alpha * t * grad_norm) {
      t <- beta * t
      par_new <- par - t * gr
    }
    par <- par_new
  }
  
  if (i == maxit)  warning("Maximal number, ", maxit, ", of iterations reached")
  
  par
}


##### Objective function #####################################

H <- function(par) {
  sum(y - f(par, x))^2
} 

grad_gd <- function(par) 1 / length(x) * sum(grad(par, 1:length(x)))    


##### Tracer #####################################

SGD_tracer <- tracer(c("par", "k"), Delta = 0) 

squared_error_mult <- function(alpha, beta, gamma, rho){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) H(par))
}

SGD_trace <- function(trace) {
  transform(
    trace,
    loss = squared_error_mult(par.1, par.2, par.3, par.4) #,
    #H_distance = abs(H(test, parameters[1], parameters[2], parameters[3], parameters[4]) - 
    #                   squared_error(test, par.1, par.2, par.3, par.4))
    )
}

###### S3 Classes ###########################################

# Parameters class

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

sim <- function(object, N, omega = 1, grid = FALSE) {
  if(grid){
    grid_sample(object$par, N)
  }
  gauss_sample(N, object$par, omega)
}

# SGD class

SGD <- function(par0, grad, n, gamma, maxiter = 100, 
                sampler = sample, cb = SGD_tracer$tracer) {
  structure(
    list(
      est = sgd(par0, grad, n, gamma, maxiter, sampler, cb),
      trace = summary(SGD_tracer),
      start_par = par0),
    class = "My_SGD"
  )
}

summary.My_SGD <- function(object) {
  SGD_trace(object$trace)
}

plot.My_SGD <- function(x) {}


