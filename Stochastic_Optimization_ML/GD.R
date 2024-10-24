## Gradient Descent #####################################

grad_desc <- function(
    par,
    grad,
    H,
    t0 = 5e-2,
    maxit = 1200,
    cb = NULL,
    epsilon = 1e-3,
    beta = 0.8,
    alpha = 0.1,
    x,
    y,
    ...) {
  for (i in 1:maxit) {
    
    n <- length(x)
    
    # Calculations of objective and gradient
    value <- H(par, x, y)
    gr <- 1 / n * grad(par, x, y)
    
    grad_norm <- sum(gr^2)
    
    # Callback
    if (!is.null(cb)) cb()
    
    t <- t0
    # Proposed descent step
    par_new <- par - t * gr
    
    # Convergence criterion 
    if (all(abs(par_new - par) <= epsilon)) break
    
    # Backtracking line search
    while (H(par_new, x, y) > value - alpha * t * grad_norm) {
      t <- beta * t
      par_new <- par - t * gr
    }
    par <- par_new
  }
  
  if (i == maxit)  warning("Maximal number, ", maxit, ", of iterations reached")
  
  par
}

grad_gd <- function(par, x, y) {
  rowSums(sapply(seq_along(x), function(i) grad(par, x[i], y[i])))
}

##### Tracer #####################################

GD_tracer <- tracer(c("par", "i"), Delta = 0) 

###### GD class ###########################################

GD <- function(par,
               grad = grad_gd,
               H = H,
               t0 = 5e-2,
               maxit = 1200,
               cb = GD_tracer$tracer,
               epsilon = 1e-3,
               beta = 0.8,
               alpha = 0.1,
               true_par = NULL,
               ...) {
  output = structure(
    list(
      est = grad_desc(par = par, 
                      grad = grad, 
                      t0 = t0,
                      maxit = maxit,
                      cb = cb,
                      epsilon = epsilon,
                      H = H,
                      beta = beta,
                      alpha = alpha,
                      ...),
      trace = summary(GD_tracer),
      start_par = par,
      true_par = true_par,
      additional_args = list(...)),
    class = "My_GD"
  )
  GD_tracer$clear()
  return(output)
}

# Summary method
summary.My_GD <- function(object) {
  object$trace
}


# Print method
print.My_GD <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("True parameters:\n")
  print(object$true_par)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,5])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,6])
}

# Plot method
plot.My_GD <- function(object, plot_no = 1, ...) {
  x <- object$additional_args$x
  y <- object$additional_args$y
  true_par <- object$true_par
  
  loss <- H_mult(x = x, y = y, 
                             alpha = object$trace$par.1, 
                             beta = object$trace$par.2,
                             gamma = object$trace$par.3,
                             rho = object$trace$par.4)
  true_par <- object$additional_args$true_par
  H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
  abs_dist_from_par <- apply(object$trace[,1:4], 1, FUN = function(par_est) 
    sum(abs(par_est - true_par)))
  
  GD_plot_df <- data.frame(object$trace, loss, H_distance)
  
  if (plot_no == 1) {
    ggplot(GD_plot_df, aes(x = .time, y = loss)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Loss vs Time", x = "Time", y = "Loss")
  } else if (plot_no == 2){
    ggplot(GD_plot_df, aes(x = .time, y = H_distance)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Distance to True Loss vs Time", x = "Time", y = "Distance") 
  } else if (plot_no == 3){
    ggplot(GD_plot_df, aes(x = .time, y = abs_dist_from_par)) +
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


plot_data.My_GD <- function(object) {
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
  
  GD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, abs_dist_from_par)
  
  return(GD_plot_df)
}

