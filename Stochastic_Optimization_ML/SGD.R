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
    gamma, # Decay schedule or a fixed learning rate
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    epoch = NULL,
    x,
    y,
    ...) {
  
  n <- length(x)
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter) {
    if (!is.null(cb)) cb()
    samp <- sampler(n)
    
    if (is.null(epoch)){
      #par <- vanilla(par, i, samp, gamma, n, ...)
      for (j in 1:n) {
        i <- samp[j]
        par <- par - gamma[k] * grad(par, x[i], y[i])
      }
    } else {
      par <- epoch(par, samp, gamma[k], ...)
    }
  }
  par
}


##### Tracer #####################################

SGD_tracer <- tracer(c("par", "k"), Delta = 0) 

###### SGD class ###########################################

SGD <- function(par0, grad, gamma, maxiter = 100, epoch = NULL,
                sampler = sample, cb = SGD_tracer$tracer, ...) {
  structure(
    list(
      est = sgd(par = par0, 
                grad = grad, 
                gamma = gamma, 
                maxiter = maxiter, 
                sampler = sample, 
                cb = cb,
                epoch = epoch,
                ...),
      trace = summary(SGD_tracer),
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
  
  loss <- H_mult(x = x, y = y, 
                 alpha = object$trace$par.1, 
                 beta = object$trace$par.2,
                 gamma = object$trace$par.3,
                 rho = object$trace$par.4)
  
  if ("true_par" %in% names(object$additional_args)) {
    true_par <- object$additional_args$true_par
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    abs_dist_from_par <- apply(object$trace[,1:4], 1, 
                               FUN = function(par_est) sum(abs(par_est - true_par)))
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
  
  loss <- H_mult(x = x, y = y, 
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

