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
    N, # Sample size
    gamma, # Decay schedule or a fixed learning rate
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    ...) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter) {
    if (!is.null(cb)) cb()
    samp <- sampler(N)
    for (j in 1:N) {
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
  grid <- exp(seq(1, 15, 1))
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

H <- function(x, y, par) {
  sum(y - f(x = x, par = par))^2
}

grad_gd <- function(par) 1 / length(x) * sum(grad(par, 1:length(x)))    


##### Tracer #####################################

SGD_tracer <- tracer(c("par", "n"), Delta = 0)

squared_error_mult <- function(x, y, alpha, beta, gamma, rho){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) H(x = x, y = y, par))
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
    grid_sample(N, object$par)
  }
  gauss_sample(N, object$par, omega)
}

#-------------------------------------------------------------------------------
#                                   SGD class                                  #
#-------------------------------------------------------------------------------

# SGD class
SGD <- function(par0, N, gamma, epoch = NULL, maxiter = 100,
                sampler = sample, cb = NULL,...) {
  structure(
    list(
      est = sgd(par0 = par0, N = N, gamma = gamma, epoch = epoch, maxiter = maxiter, sampler = sampler, cb = cb,...),
      trace = summary(cb),
      start_par = par0,
      additional_args = list(...)),
    class = "My_SGD"
  )
}


# Summary method
summary.My_SGD <- function(object) {
  return(object$trace)
}


# Print method
print.My_SGD <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,5])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,6])
}


# Plot method
plot.My_SGD <- function(object, plot_no = 1, ...) {
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
  
  SGD_plot_df <- data.frame(object$trace, loss, H_distance)
  
  if (plot_no == 1) {
    ggplot(SGD_plot_df, aes(x = .time, y = loss)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Loss vs Time", x = "Time", y = "Loss")
  } else if (plot_no == 2){
    ggplot(SGD_plot_df, aes(x = .time, y = H_distance)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Distance to True Loss vs Time", x = "Time", y = "Distance") 
  } else if (plot_no == 3){
    ggplot(SGD_plot_df, aes(x = .time, y = abs_dist_from_par)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Sum of absolute distance to true parameters vs Time", x = "Time", y = "Distance") 
  } else {
    stop("Invalid plot number. Please choose 1 2 or 3.")
  }
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


