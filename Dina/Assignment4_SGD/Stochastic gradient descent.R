f <- function(x, par){
  alpha <- par[[1]]
  beta <- par[[2]]
  gamma <- par[[3]]
  rho <- par[[4]]
  
  return(gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha)))
}



gradient <- function(par, x, y, ...){
  alpha <- par[[1]]
  beta <- par[[2]]
  gamma <- par[[3]]
  rho <- par[[4]]
  
  
  expbetalogxalpha <- exp(beta * log(x) - alpha)
  
  identical_part <- - 2 * (y - f(x = x, par = par))
  
  grad_alpha <- mean(identical_part * (rho - gamma) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_beta <- - mean(identical_part * (rho - gamma) * log(x) * expbetalogxalpha / (1 + expbetalogxalpha)^2)
  grad_gamma <- mean(identical_part * (1 - 1 / (1 + expbetalogxalpha)))
  grad_rho <- mean(identical_part / (1 + expbetalogxalpha))
  
  return(c(grad_alpha, grad_beta, grad_gamma, grad_rho))
}








sgd <- function(
    par0,
    x,
    y,
    grad,
    gamma, # Decay schedule or a fixed learning rate
    epoch = vanilla,
    maxiter = 100, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    epsilon = 1e-6, # Tolerance for convergence
    ...) {
  
  if (is.function(gamma)) gamma <- gamma(1:maxiter)
  gamma <- rep_len(gamma, maxiter)
  
  N <- length(x)
  par <- par0
  
  for (n in 1:maxiter) {
    
    if (!is.null(cb)) cb$tracer()
    samp <- sampler(N)
    par <- epoch(par = par, samp = samp, gamma = gamma[n], 
                 grad = grad, x = x, y = y, ...)
    if (sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon))
      break
  }
  return(par)
}



##### Objective function #####################################

H <- function(x, y, par) {
  mean((y - f(x = x, par = par))^2)
}

H_mult <- function(x, y, param){
  apply(param, 1, function(par) H(x = x, y = y, par))
}




###### S3 Classes ###########################################

#-------------------------------------------------------------------------------
#                                   SGD class                                  #
#-------------------------------------------------------------------------------

SGD_tracer <- tracer(c("par", "n"), Delta = 0)

# SGD class
SGD <- function(par0, x, y, grad, gamma, epoch = vanilla, maxiter = 100,
                sampler = sample, cb = SGD_tracer, epsilon = 1e-6, true_par = NULL, ...) {
  
  
  if(!is.null(cb)) cb$clear()
  
  output <- structure(
    list(x = x, y = y,
      est = sgd(par0 = par0,
                x = x,
                y = y,
                grad = grad,
                gamma = gamma, # Decay schedule or a fixed learning rate
                epoch = epoch,
                maxiter = maxiter, # Max epoch iterations
                sampler = sampler, # How data is resampled. Default is a random permutation
                cb = cb,
                epsilon = epsilon, # Tolerance for convergence
                ...),
      trace = summary(cb),
      start_par = par0,
      epoch = epoch,
      true_par = true_par,
      gamma = gamma,
      additional_args = list(...)),
    class = "My_SGD"
  )
  
  if(!is.null(cb)) cb$clear()
  
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
#test <- SGD(par0 = param$par, grad = gradient, gamma = 0.01, x = x, y = y, true_par = param$par)
#test$trace

# Plot method
plot.My_SGD <- function(object, plot_no = 1, ...) {
  x <- object$x
  y <- object$y
  
  loss <- H_mult(x = x, y = y, param = object$trace[,1:4])
  true_par <- object$true_par
  
  # Requires true_par to be specified
  if (!is.null(true_par)){
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    DistToTruepar <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(sqrt((par_est - true_par)^2)))
    
    SGD_plot_df <- data.frame(object$trace, loss, H_distance, DistToTruepar)
  } else {
    SGD_plot_df <- data.frame(object$trace, loss)
  }
  
  
  
  if (plot_no == 1) {
    ggplot(SGD_plot_df, aes(x = .time, y = loss)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Loss vs Time", x = "Time", y = "Loss")
  } else if (plot_no == 2){
    
    if (is.null(true_par)) stop("True parameters must be specified as true_par")
    
    ggplot(SGD_plot_df, aes(x = .time, y = H_distance)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Suboptimality", x = "Time", y = expression("|H("*theta^"'"*") - H("*theta*")|")) 
    
  } else if (plot_no == 3){
    
    if (is.null(true_par)) stop("True parameters must be specified as true_par")
    
    ggplot(SGD_plot_df, aes(x = .time, y = DistToTruepar)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Distance to true parameters", x = "Time", y = "Distance") 
    
  } else {
    stop("Invalid plot number. Please choose 1 2 or 3.")
  }
}


# Method to extract plot data
plot_data <- function(x) {
  UseMethod("plot_data")
}


plot_data.My_SGD <- function(object) {
  x <- object$x
  y <- object$y
  true_par <- object$true_par
  
  loss <- H_mult(x = x, y = y, object$trace[,1:4])
  
  if (!is.null(true_par)){
    true_par <- object$true_par
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    DistToTruepar <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(sqrt((par_est - true_par)^2)))
    
    SGD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, DistToTruepar)
  } else {
    SGD_plot_df <- data.frame(".time" = object$trace$.time, loss)
  }
  
  
  return(SGD_plot_df)
}



