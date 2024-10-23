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

## Function  ###########################################

f <- function(par, x) {
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  gamma + (rho - gamma)/(1 + exp(beta * log(x) - alpha))
}


## Gradient  ###########################################


nabla_f <- function(par, x){
    
    alpha <- par[1]
    beta <- par[2]
    gamma <- par[3]
    rho <- par[4]
    
    const_term <- exp(beta * log(x) - alpha)
    
    d_alpha <- (rho - gamma) / (1 + const_term)^2 * const_term
    d_beta <- - log(x) * (rho - gamma) / (1 + const_term)^2 * const_term
    d_gamma <- 1 - 1 / (1 + const_term)
    d_rho <- 1 / (1 + const_term)
    
    return(c(d_alpha, d_beta, d_gamma, d_rho))
}

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

##### Objective function #####################################

H <- function(par, x, y) {
  mean((y - f(par, x))^2)
} 

H_mult <- function(alpha, beta, gamma, rho, x, y){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) H(par, x, y))
}

## Sampler  ##############################################

gauss_sample <- function(N, par, omega = 1) {
  log_x <- rnorm(N, 0, omega^2)
  x <- exp(log_x)
  y <- f(par, x) + rnorm(N, 0, 0.5)
  return(data.frame(x = x, y = y))
}

grid_sample <- function(N, par) {
  grid <- exp(1:15)
  x <- sample(grid, N, replace = TRUE)
  y <- f(par, x) + rnorm(N, 0, 0.5)
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
    data <- grid_sample(N, object$par)
  }
  data <- gauss_sample(N, object$par, omega)
  #For scaling
  if(scale){
    scaled_data <- scale(data) 
    # Shift the data so that the minimum value is 1 (or any positive value)
    data <- scaled_data + abs(min(scaled_data)) + 1
  }
  return(data)
}
