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

## Parameter Class ###############################################

generate <- function(n = 2000, par) {
  data_list <- list()
  for(i in 1:length(par$mu)){
    w <- rchisq(n, df = par$nu)
    x <- rnorm(n, mean = par$mu, sd = sqrt(par$sigma2 * par$nu / w))
    data_list[[i]] <- list(x = x, w = w)
  }
  return(data_list)
}

parameters <- function(mu, sigma2, nu) {
  data <- generate(n = 4000,
                   par = list(mu = mu, sigma2 = sigma2, nu = nu))
  structure(
    list(
      mu = mu,
      sigma2 = sigma2,
      nu = nu, 
      data = data),
    class = "My_parameters"
    )
}

# Plot function

plot.My_parameters <- function(x, ...) {
  data <- x$data
  n <- length(data)
  par(mfrow = c(1, n))
  for(i in 1:n){
    hist(data[[i]]$x, main = paste("Component", i), xlab = "x")
  }
}


## Complete data MLE estimators ##################################################

mle.mu <- function(x, w){
  sum(x*w)/(sum(w))
}

mle.sigma2 <- function(x, w, mu_mle, nu){
  sum(w*(x-mu_mle)^2) / (nu * length(x))
}

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
