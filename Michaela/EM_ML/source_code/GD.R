# Gradient of negative log likelihood
grad <- function(mu, sigma, nu, x){
  n <- length(x)
  d_sigma <- n/sigma - (nu + 1) / sigma * sum((x - mu)^2 / (nu * sigma^2 + (x - mu)^2))
  d_mu <- - (nu + 1) * sum((x - mu) / (nu * sigma^2 + (x - mu)^2))
  return(c(d_mu, d_sigma))
}

# Negative loglikelihood
neg_loglik <- function(par, x, nu) {
  - loglik(par = par, x = x, nu)
}


# Gradient descent algorithm 
grad_desc <- function(
    par,
    grad,
    H,
    t0 = 1,
    maxit = 1200,
    cb = NULL,
    epsilon = 1e-6,
    beta = 0.8,
    alpha = 0.1,
    x,
    nu,
    ...) {
  
  n <- length(x)
  
  for (i in 1:maxit) {
    
    #browser()
    
    # Calculations of objective and gradient
    value <- 1 / n * H(x = x, par = par, nu = nu)
    gr <- 1/n * grad(mu = par[1], sigma = sqrt(par[2]), nu = nu, x = x)
    
    grad_norm <- sum(gr^2)
    
    t <- t0
    # Proposed descent step
    par_new <- par - t * gr
    
    # Backtracking line search
    while (1 / n * H(x = x, par = par_new, nu = nu) > value - alpha * t * grad_norm) {
      t <- beta * t
      par_new <- par - t * gr
    }
    
    # Callback
    if (!is.null(cb)) cb()
    
    # Convergence criterion 
    if (sum((par_new - par)^2) <= epsilon * (sum(par_new^2) + epsilon)) break
    
    par <- par_new
  }
  
  if (i == maxit)  warning("Maximal number, ", maxit, ", of iterations reached")
  
  par
}

##### Tracer #####################################

GD_tracer <- tracer(c("par_new", "par", "value", "gr", "grad_norm", "i"), Delta = 0) 


###### GD class ###########################################

GD <- function(par,
               grad = grad,
               H = neg_loglik,
               t0 = 1,
               maxit = 1200,
               cb = GD_tracer,
               epsilon = 1e-6,
               beta = 0.8,
               alpha = 0.1,
               nu = NULL,
               ...) {
  output = structure(
    list(
      est = grad_desc(par = par, 
                      grad = grad, 
                      t0 = t0,
                      maxit = maxit,
                      cb = cb$tracer,
                      epsilon = epsilon,
                      H = H,
                      beta = beta,
                      alpha = alpha,
                      nu = nu,
                      ...),
      trace = summary(GD_tracer),
      start_par = par,
      nu = nu,
      additional_args = list(...)),
    class = "My_GD"
  )
  GD_tracer$clear()
  return(output)
}


# Print method
print.My_GD <- function(object){
  cat("Estimated parameters:\n")
  print(object$est)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,9])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,10])
}


# Plot method
plot.My_GD <- function(object, plot_no = 1, ...) {
  x <- object$x
  nu <- object$nu
  
  GD_plot_df <- data.frame(object$trace)
  
  if (plot_no == 1) {
    ggplot() +
      geom_line(data = GD_plot_df, aes(x = .time, y = value, color = "GD")) +
      scale_y_log10() +
      labs(title = "Neg. Loglik vs Time", x = "Time", y = "Neg. Loglik")
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

