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
    maxiter = 150, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    epoch = vanilla,
    m = 1, # Batch size
    x,
    y,
    ...) {
  
  n <- length(x)
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  
  for (k in 1:maxiter) {
    
    if (!is.null(cb)) cb()
    samp <- sampler(n)
    par <- epoch(par = par, samp = samp, gamma = gamma[k], 
                 grad = grad, n = n, x = x, y = y, m = m)
    
  }
  par
}

##### Tracer #####################################

SGD_tracer <- tracer(c("par", "k"), Delta = 0) 

###### SGD class ###########################################

SGD <- function(par0, grad, gamma, maxiter = 150, epoch = vanilla,
                sampler = sample, cb = SGD_tracer, m = 1, 
                true_par = NULL, ...) {
  output = structure(
    list(
      est = sgd(par = par0, 
                grad = grad, 
                gamma = gamma, 
                maxiter = maxiter, 
                sampler = sample, 
                cb = cb$tracer,
                epoch = epoch,
                m = m,
                ...),
      trace = summary(cb),
      start_par = par0,
      true_par = true_par,
      additional_args = list(...)),
    class = "My_SGD"
  )
  cb$clear()
  return(output)
}

# Summary method
summary.My_SGD <- function(object) {
  return(object$trace)
}

# Print method
print.My_SGD <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("True parametes:\n")
  print(object$true_par)
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
  
  true_par <- object$true_par
  H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
  abs_dist_from_par <- apply(object$trace[,1:4], 1, 
                               FUN = function(par_est) sum(abs(par_est - true_par)))
  
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
  
  true_par <- object$true_par
  H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
  abs_dist_from_par <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(abs(par_est - true_par)))
  
  SGD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, abs_dist_from_par)
  
  return(SGD_plot_df)
}

##### Decay Scheduler #####################################

vanilla <- function(par, samp, gamma, grad, n, x, y, ...){
  for (j in 1:n) {
    i <- samp[j]
    par <- par - gamma * grad(par, x[i], y[i])
  }
  return(par)
}

batch <- function(par, samp, gamma, grad, n, x, y, m, ...){
  M <- floor(n / m)
  for (j in 0:(M - 1)) {
    i <- samp[(j * m + 1):(j * m + m)]
    par <- par - 1 / m * gamma * grad(par, x[i], y[i])
  }
  return(par)
}

decay_scheduler <- function(gamma0 = 1, # Initial learning rate
                            a = 1, 
                            K = 1, 
                            gamma1,     # Optional target learning rate
                            n1         # We want to reach after n1 iterations
                            ){ 
  force(a)
  if (!missing(gamma1) && !missing(n1))
    K <- n1^a * gamma1 / (gamma0 - gamma1)
  b <- gamma0 * K
  function(n) b / (K + n^a)
}


adam <- function() {
  rho <- v <- 0
  function(
   #    par,
    samp,
    gamma,
    grad,
    m = 50,         # Mini-batch size
    beta1 = 0.9,    # Momentum memory
    beta2 = 0.9,    # Momentum memory
    ...
    
  ){
    M <- floor(length(samp) / m) 
    for (j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      gr <- grad(par, i, ...)
      rho <<- beta1 * rho + (1 - beta1) * gr
      v <<- beta2 * v + (1 - beta2) * gr^2
      par <- par - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par
  } 
}


momentum <- function() {
  rho <- 0
  function(
    par,
    samp,
    gamma,
    grad,
    m = 50,         # Mini-batch size
    beta = 0.95,    # Momentum memory
    ...
  ){
    M <- floor(length(samp) / m) 
    for (j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      # Using '<<-' assigns the value to rho in the enclosing
      environment
      rho <<- beta * rho + (1 - beta) * grad(par, i, ...)
      par <- par - gamma * rho
    }
    par
  } 
}

rate_momentum <- decay_scheduler(gamma0 = 1, a = 1, gamma1 = 1e-1)
rate_adam <- decay_scheduler(gamma0 = 1e-1, a = 1, gamma1 = 1e-5)


########### Rcpp ############################

library(Rcpp)

cppFunction('
NumericVector gradient_rcpp(NumericVector par, NumericVector indices, NumericVector x, NumericVector y) {
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double gamma = par[2];
  double rho = par[3];

  // Initialize gradients
  double grad_alpha = 0.0;
  double grad_beta = 0.0;
  double grad_gamma = 0.0;
  double grad_rho = 0.0;

  int n = indices.size(); // Number of indices

  // Loop over the indices
  for (int idx = 0; idx < n; ++idx) {
    int i = indices[idx] - 1;  // Convert from r to c++ indexing

    // Get individual data point
    double x_i = x[i];
    double y_i = y[i];

    // Calculating f(x_i, par)
    double f_x_i = gamma + (rho - gamma) / (1 + exp(beta * log(x_i) - alpha));

    // Exponential term
    double expbetalogxalpha = exp(beta * log(x_i) - alpha);

    // Identical part used in gradients
    double identical_part = -2 * (y_i - f_x_i);

    // Accumulate gradients for all indices
    grad_alpha += (identical_part * (rho - gamma) * expbetalogxalpha) 
    / pow(1 + expbetalogxalpha, 2);
    grad_beta += -(identical_part * (rho - gamma) * log(x_i) * expbetalogxalpha) 
    / pow(1 + expbetalogxalpha, 2);
    grad_gamma += identical_part * (1 - 1 / (1 + expbetalogxalpha));
    grad_rho += identical_part / (1 + expbetalogxalpha);
  }

  // Return the mean of accumulated gradients
  return NumericVector::create(grad_alpha / n, grad_beta / n, grad_gamma / n, grad_rho / n);
}
')
