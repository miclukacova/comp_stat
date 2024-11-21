# Gradient clipping
grad_clipped <- function(grad, grad_norm, threshold = 1){
  if (grad_norm > threshold){
    grad <- grad / grad_norm * 2 * threshold
  }
  return(grad)
}





# Gradient descent
grad_desc <- function(
    par0,
    H,
    grad,
    d = 0.8,
    c = 0.1,
    gamma0 = 1,
    epsilon = 1e-9,
    maxiter = 1000,
    beta = 0.9,
    cb = NULL,
    clipping = FALSE,
    momentum = FALSE,
    threshold = 1,
    ...) {
  
    gr <- 0
    for (i in 1:maxiter) {
      
      
      
      # Calculations of objective and gradient
      value <- H(par0)
      if (momentum){
        gr <- beta * gr + (1 - beta) * grad(par0)
      } else{
        gr <- grad(par0)
      }
      
      
      # Gradient clipping
      if (clipping) {
        grad_norm <- sum(gr^2)
        if (grad_norm > threshold){
          gr <- gr / grad_norm * threshold
        }
      }
      
      grad_norm <- sum(gr^2)
      
      # Convergence criterion based on gradient norm
      if (grad_norm <= epsilon) break
      
      
      gamma <- gamma0
      # Proposed descent step
      par1 <- par0 - gamma * gr
      
      
      # Backtracking line search
      while (H(par1) > value - c * gamma * grad_norm) {
        gamma <- d * gamma
        par1 <- par0 - gamma * gr
      }
      
      # Callback
      if(!is.null(cb)) cb$tracer()
      
      par0 <- par1
    }
    
    if (i == maxiter)  warning("Maximal number, ", maxiter, ", of iterations reached")
    
    return(par0)
}




# Objects
GD_tracer <- tracer(c("par0", "par1", "value", "gr", "grad_norm", "i"), Delta = 0)

# EM class
GD <- function(x, y,
               par0,
               H,
               grad,
               d = 0.8,
               c = 0.1,
               gamma0 = 0.1,
               epsilon = 1e-9,
               maxiter = 1000,
               beta = 0.9,
               cb = GD_tracer,
               threshold = 1,
               true_par = NULL,
               momentum = FALSE,
               clipping = FALSE,...) {
  
  if(!is.null(cb)) cb$clear()
  
  
  # Function factory
  GD_function_factory_cpp <- function(x, y){
    
    # Objective function
    H_par <- function(par){
      H(x = x, y = y, par = par)
    }
    
    # Gradient
    grad_par <- function(par){
      grad(x = x, y = y, par = par)
    }
    
    return(list(H = H_par, grad = grad_par))
  }
  
  GD_func_fac <- GD_function_factory_cpp(x = x, y = y)
  
  H_mult <- function(par){
    H <- GD_func_fac$H
    return(apply(par, 1, H))
  }
  
  
  output <- structure(
    list(
      x = x,
      y = y,
      est = grad_desc(par0 = par0,
                      H = GD_func_fac$H,
                      grad = GD_func_fac$grad,
                      d = d,
                      c = c,
                      gamma0 = gamma0,
                      epsilon = epsilon,
                      maxiter = maxiter,
                      beta = beta,
                      momentum = momentum,
                      cb = cb,
                      threshold = threshold,
                      clipping = clipping,...),
      trace = summary(cb),
      gamma0 = gamma0,
      par0 = par0,
      H = GD_func_fac$H,
      grad = GD_func_fac$grad,
      H_mult = H_mult,
      true_par = true_par,
      clipping = clipping,
      momentum = momentum,
      beta = beta,
      d = d,
      c = c,
      additional_args = list(...)),
    class = "My_GD"
  )
  
  if(!is.null(cb)) cb$clear()
  
  return(output)
}



print.My_GD <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("Number of data points:\n")
  print(length(object$x))
  cat("Initial stepsize:\n")
  print(object$gamma0)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,15])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,16])
  cat("Momentum:\n")
  print(object$momentum)
  cat("Gradient Clipping:\n")
  print(object$clipping)
}


plot.My_GD <- function(object, plot_no = 1, ...) {
  x <- object$x
  y <- object$y
  
  loss <- H_mult(x = x, y = y, param = object$trace[,1:4])
  true_par <- object$true_par
  
  # Requires true_par to be specified
  if (!is.null(true_par)){
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    DistToTruepar <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(sqrt((par_est - true_par)^2)))
    
    GD_plot_df <- data.frame(object$trace, loss, H_distance, DistToTruepar)
  } else {
    GD_plot_df <- data.frame(object$trace, loss)
  }
  
  
  
  if (plot_no == 1) {
    ggplot(GD_plot_df, aes(x = .time, y = loss)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Loss vs Time", x = "Time", y = "Loss")
  } else if (plot_no == 2){
    
    if (is.null(true_par)) stop("True parameters must be specified as true_par")
    
    ggplot(GD_plot_df, aes(x = .time, y = H_distance)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Suboptimality", x = "Time", y = expression("|H("*theta^"'"*") - H("*theta*")|")) 
    
  } else if (plot_no == 3){
    
    if (is.null(true_par)) stop("True parameters must be specified as true_par")
    
    ggplot(GD_plot_df, aes(x = .time, y = DistToTruepar)) +
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


plot_data.My_GD <- function(object) {
  x <- object$x
  y <- object$y
  true_par <- object$true_par
  
  loss <- H_mult(x = x, y = y, object$trace[,1:4])
  
  if (!is.null(true_par)){
    true_par <- object$true_par
    H_distance <- abs(H(x = x, y = y, par = true_par) - loss)
    DistToTruepar <- apply(object$trace[,1:4], 1, FUN = function(par_est) sum(sqrt((par_est - true_par)^2)))
    
    GD_plot_df <- data.frame(".time" = object$trace$.time, loss, H_distance, DistToTruepar)
  } else {
    GD_plot_df <- data.frame(".time" = object$trace$.time, loss)
  }
  
  
  return(GD_plot_df)
}



