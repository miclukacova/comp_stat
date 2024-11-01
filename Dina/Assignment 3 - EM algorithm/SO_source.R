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

## EM ###############################################

EM_alg <- function(x, param, max_iter = 20, epsilon = 1e-10, cb = NULL){
  mu_mark <- param$mu #Getting initial parameters
  sigma_mark <- param$sigma
  ny <- param$ny
  k <- (ny + 1) / 2   #Defining i and k, which only depends on ny:
  i <- 0
  
  while (i < max_iter) {
    mu_old <- mu_mark
    sigma_old <- sigma_mark
    t_old <- 2 / (1 + (x - mu_old)^2 / (ny * sigma_old))
    
    mu_mark <- sum(t_old * x) / sum(t_old)
    sigma_mark <- 1 / (length(x) * ny) * k * sum(t_old * (x - mu_mark)^2)
    #Calling  cb
    if(!is.null(cb)) cb()
    #Checking convergence
    if (sum((c(mu_mark, sigma_mark) - c(mu_old, sigma_old))^2) 
        <= epsilon * (sum(c(mu_mark, sigma_mark)^2) + epsilon)) break
    i <- i + 1
  }
  c(mu_mark, sigma_mark)
}

## log-likelihood  ###########################################

log_lik <- function(m, s, ny, x = samp$x){
  k <- (ny + 1) / 2
  return( - length(x) * log(sqrt(s)) - sum(k * log(1 + (x - m)^2 / (ny * s))))
}

## Norm of the difference between parameters ###########################################

par_norm_diff <- function(mu_old, sigma_old, mu_mark, sigma_mark){
  return(sqrt((mu_old - mu_mark)^2 + (sigma_old - sigma_mark)^2))
}

## Gradient  ###########################################

grad <- function(m, s, ny, x){
  k <- (ny + 1) / 2
  t <- sum(2 / (1 + (x - m)^2 / (ny * s)))
  grad_m <- - t * sum(k * (x - m) / (ny * s))
  grad_s <- - sum(1 / s) + t * sum(k * (x - m)^2 / (ny * sqrt(s)^3))
  return(c(grad_m, grad_s))
}

## Sampler  ##############################################

samp_our_distribution <- function(n, param){
    w <- rchisq(n, df = param$ny)
    x <- rnorm(n, mean = param$mu, sd = sqrt((param$ny * param$sigma) / w))
    y <- cbind(x, w)
    return(y)
  }

## "Theoretical parameters" #####################################

theo_par <- function(x, w, param){
  mu_opt <- sum(w * x) / sum(w)
  sigma_opt <- 1/(length(x) * param$ny) * sum(w * (x - mu_opt)^2)
  return(data.frame(mu = mu_opt, sigma = sigma_opt ))
}

##### Fischer Information #####################################



##### Tracer #####################################

EM_tracer <- tracer(c("mu_mark", "sigma_mark", "mu_old", "sigma_old", "ny", "i"), Delta = 0)

log_lik_mult <- function(m, s, ny, x){
  param <- cbind(m, s, ny)
  print(param)
  apply(param, 1, function(par) log_lik(par[1], par[2], par[3], x = x))
}

grad_mult <- function(m, s, ny, x){
  param <- cbind(m, s, ny)
  apply(param, 1, function(par) grad(par[1], par[2], par[3], x = x))
}

###### S3 Classes ###########################################

# Parameters class

parameters <- function(mu, sigma, ny) {
  structure(
    list(
      mu = mu,
      sigma = sigma,
      ny = ny, 
      par = c(mu, sigma, ny)),
    class = "My_params"
  )
}

sim <- function(x) {
  UseMethod("sim")
}

sim <- function(object, N) {
  return(as.data.frame(samp_our_distribution(n = N, param = object)))
}


#-------------------------------------------------------------------------------
#                                   EM class                                  #
#-------------------------------------------------------------------------------

# EM class
EM <- function(x, par0, max_iter, epsilon, cb = NULL,...) {
  structure(
    list(
      est = EM_alg(x = x, param = par0, max_iter = max_iter, epsilon = epsilon, cb = cb$tracer,...),
      trace = summary(cb),
      start_par = par0,
      x = x,
      additional_args = list(...)),
    class = "My_EM"
  )
}


# Summary method
summary.My_EM <- function(object) {
  return(object$trace[c(1,2,5,6,7)])
}


# Print method
print.My_EM <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,6])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,7])
}

#Plot method for convergence plots
plot.My_EM <- function(object, time = F, ...) {
  x <- object$x
  
  par_diff <- par_norm_diff(object$trace$mu_old, object$trace$sigma_old, object$trace$mu_mark, object$trace$sigma_mark)
  
  ll <- log_lik_mult(object$trace$mu_mark, object$trace$sigma_mark, object$trace$ny, x)
  log_lik_convergence <- 
}


# Method to extract plot data
plot_data <- function(x) {
  UseMethod("plot_data")
}


plot_data.My_SGD <- function(object) {
  x <- object$additional_args$x
  y <- object$additional_args$y
  
  loss <- squared_error_mult(x = x, y = y, 
                             alpha = object$trace$par.1, 
                             beta = object$trace$par.2,
                             gamma = object$trace$par.3,
                             rho = object$trace$par.4)
  
  if ("true_par" %in% names(object$additional_args)) {
    true_par <- object$additional_args$true_par
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    abs_dist_from_par <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(abs(par_est - true_par)))
  }
  
  SGD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, abs_dist_from_par)
  
  return(SGD_plot_df)
}


