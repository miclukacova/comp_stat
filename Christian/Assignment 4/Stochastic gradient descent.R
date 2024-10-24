f <- function(x, par){
  alpha <- par[[1]]
  beta <- par[[2]]
  gamma <- par[[3]]
  rho <- par[[4]]
  
  return(gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha)))
}



gradient <- function(par, i, x, y,...){
  alpha <- par[[1]]
  beta <- par[[2]]
  gamma <- par[[3]]
  rho <- par[[4]]
  
  
  x_i <- x[i]
  y_i <- y[i]
  
  expbetalogxalpha <- exp(beta * log(x_i) - alpha)
  
  identical_part <- - 2 * (y_i - f(x = x_i, par = par))
  
  grad_alpha <- mean(identical_part * (rho - gamma) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_beta <- - mean(identical_part * (rho - gamma) * log(x_i) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_gamma <- mean(identical_part * (1 - 1 / (1 + expbetalogxalpha)))
  grad_rho <- mean(identical_part / (1 + expbetalogxalpha))
  
  return(c(grad_alpha, grad_beta, grad_gamma, grad_rho))
}




sgd <- function(
    par0,
    N, # Sample size
    gamma, # Decay schedule or a fixed learning rate
    epoch = NULL,
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    ...) {
  
  if (is.function(gamma)) gamma <- gamma(1:maxiter)
  gamma <- rep_len(gamma, maxiter)
  
  par <- par0
  
  for (n in 1:maxiter) {
    if (!is.null(cb)) cb$tracer()
    
    samp <- sampler(N)
    
    if (is.null(epoch)){
      par <- vanilla(par, i, samp, gamma[n], n = n, ...)
    } else {
      par <- epoch(par, samp, gamma[n], ...)
    }
  }
  par
}



##### Objective function #####################################

H <- function(x, y, par) {
  mean((y - f(x = x, par = par))^2)
}

squared_error_mult <- function(x, y, param){
  apply(param, 1, function(par) H(x = x, y = y, par))
}

###### S3 Classes ###########################################

#-------------------------------------------------------------------------------
#                                   SGD class                                  #
#-------------------------------------------------------------------------------

# SGD class
SGD <- function(par0, N, gamma, epoch = NULL, maxiter = 100,
                sampler = sample, cb = NULL,...) {
  output <- structure(
    list(
      est = sgd(par0 = par0, N = N, gamma = gamma, epoch = epoch, maxiter = maxiter, sampler = sampler, cb = cb,...),
      trace = summary(cb),
      start_par = par0,
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
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,5])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,6])
}


# Plot method
plot.My_SGD <- function(object, plot_no = 1, ...) {
  x <- object$additional_args$x
  y <- object$additional_args$y
  
  loss <- squared_error_mult(x = x, y = y, param = object$trace[,1:4])
  
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
  
  loss <- squared_error_mult(x = x, y = y, object$trace[,1:4])
  
  if ("true_par" %in% names(object$additional_args)) {
    true_par <- object$additional_args$true_par
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    abs_dist_from_par <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(abs(par_est - true_par)))
  }
  
  SGD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, abs_dist_from_par)
  
  return(SGD_plot_df)
}



