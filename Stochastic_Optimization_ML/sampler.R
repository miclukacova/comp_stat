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

grad <- function(par, x, y) {
  
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  d_alpha <- (rho - gamma)/(1 + exp(beta * log(x) - alpha))^2 * exp(beta * log(x) - alpha)
  d_beta <- - log(x) * (rho - gamma)/(1 + exp(beta * log(x) - alpha))^2 *
    exp(beta * log(x) - alpha)
  d_gamma <- 1 - 1/(1 + exp(beta * log(x) - alpha))
  d_rho <- 1/(1 + exp(beta * log(x) - alpha))
  
  return(- c(d_alpha, d_beta, d_gamma, d_rho) * (y - f(par, x)))
}

##### Objective function #####################################

H <- function(par, x, y) {
  sum(y - f(par, x))^2
} 

H_mult <- function(alpha, beta, gamma, rho, x, y){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) H(par, x, y))
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

sim <- function(object, N, omega = 1, grid = FALSE) {
  if(grid){
    grid_sample(object$par, N)
  }
  gauss_sample(N, object$par, omega)
}
