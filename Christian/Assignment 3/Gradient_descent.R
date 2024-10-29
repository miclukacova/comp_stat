# Calculating the marginal log likelihood

marginal_log_likelihood <- function(x, mu, sigma, nu){
  n <- length(x)
  return(-n*log(sigma) - (nu + 1)/2 * sum(log(1 + (x - mu)^2/(nu * sigma^2))))
}


# Calculating the gradient of the marginal log likelihood

d_marginal_log_likelihood <- function(x, mu, sigma, nu){
  n <- length(x)
  
  d_sigma <- -n/sigma + (nu + 1) / sigma * sum((x - mu)^2 / (nu * sigma^2 + (x - mu)^2))
  d_mu <- (nu + 1) * sum((x - mu) / (nu * sigma^2 + (x - mu)^2))
  
  return(c(d_mu, d_sigma))
}

d_marginal_log_likelihood2 <- function(x, mu, sigma, nu){
  n <- length(x)
  
  x_minus_mu <- x - mu
  x_minus_mu_sq <- x_minus_mu^2
  
  d_sigma <- -n/sigma + (nu + 1) / nu * sum(x_minus_mu_sq / (sigma^3 * (1 + x_minus_mu_sq/(nu * sigma^2))))
  d_mu <- (nu + 1) / nu * sum(x_minus_mu / (sigma^2 + x_minus_mu_sq / nu))
  
  return(c(d_mu, d_sigma))
}


d_marginal_log_likelihood2_par <- function(data, par, nu){
  n <- length(data)
  
  mu <- par[1]
  sigma <- par[2]
  
  x_minus_mu <- data - mu
  x_minus_mu_sq <- x_minus_mu^2
  
  d_sigma <- -n/sigma + (nu + 1) / nu * sum(x_minus_mu_sq / (sigma^3 * (1 + x_minus_mu_sq/(nu * sigma^2))))
  d_mu <- (nu + 1) / nu * sum(x_minus_mu / (sigma^2 + x_minus_mu_sq / nu))
  
  return(c(d_mu, d_sigma))
}



# Gradient clipping
grad_clipped <- function(grad, grad_norm, threshold = 1){
  if (grad_norm > threshold){
    grad <- grad / grad_norm * threshold
  }
  return(grad)
}

# Function factory

GD_function_factory <- function(x, nu){
  
  n <- length(x)
  
  # Negative marginal log-likelihood
  H <- function(par){
    - marginal_log_likelihood(x = x, par[[1]], par[[2]], nu = nu) / n
  }
  
  # Negative gradient of marginal log-likelihood
  grad <- function(par){
    - d_marginal_log_likelihood(x = x, par[[1]], par[[2]], nu = nu) / n
  }
  
  return(list(H = H, grad = grad))
}



GD_function_factory_cpp <- function(x, nu){
  
  n <- length(x)
  
  # Negative marginal log-likelihood
  H <- function(par){
    - marginal_log_likelihood_cpp(x = x, par[[1]], par[[2]], nu = nu) / n
  }
  
  # Negative gradient of marginal log-likelihood
  grad <- function(par){
    - d_marginal_log_likelihood_cpp(x = x, par[[1]], par[[2]], nu = nu) / n
  }
  
  return(list(H = H, grad = grad))
}



# Negative marginal log_likelihood for a vector of parameters

H_mult <- function(x, par, nu){
  H <- GD_function_factory(x, nu)$H
  return(apply(par, 1, H))
}




# Gradient descent
grad_desc <- function(
    par0,
    H,
    grad,
    d = 0.8,
    c = 0.1,
    gamma0 = 0.1,
    epsilon = 1e-4,
    maxiter = 1000,
    cb = NULL,
    clipping = FALSE,
    ...) {
  for (i in 1:maxiter) {
    
    n <- length(x)
    
    # Calculations of objective and gradient
    value <- H(par0)
    gr <- grad(par0)
    grad_norm <- sum(gr^2)
    
    if (clipping) {
      gr <- grad_clipped(grad = gr, grad_norm = grad_norm)
    }
    
    #browser()
    
    
    
    gamma <- gamma0
    # Proposed descent step
    par1 <- par0 - gamma * gr
    
    # Callback
    if(!is.null(cb)) cb$tracer()
    
    # Convergence criterion based on gradient norm
    if (grad_norm <= epsilon) break
    
    #browser()
    
    # Backtracking line search
    while (H(par1) > value - c * gamma * grad_norm) {
      gamma <- d * gamma
      par1 <- par0 - gamma * gr
    }
    par0 <- par1
  }
  
  if (i == maxiter)  warning("Maximal number, ", maxiter, ", of iterations reached")
  
  return(par1)
}







# Objects
GD_tracer <- tracer(c("par0", "par1", "value", "gr", "grad_norm", "i"), Delta = 0)

# EM class
GD <- function(x,
               nu,
               par0,
               function_factory,
               d = 0.8,
               c = 0.1,
               gamma0 = 0.1,
               epsilon = 1e-4,
               maxiter = 1000,
               cb = NULL,
               clipping = FALSE,...) {
  
  GD_func_fac <- function_factory(x, nu)
  
  
  output <- structure(
    list(
      x = x,
      nu = nu,
      est = grad_desc(par0 = par0,
                      H = GD_func_fac$H,
                      grad = GD_func_fac$grad,
                      d = d,
                      c = c,
                      gamma0 = gamma0,
                      epsilon = epsilon,
                      maxiter = maxiter,
                      cb = cb,
                      clipping = clipping,...),
      trace = summary(cb),
      par0 = par0,
      additional_args = list(...)),
    class = "My_GD"
  )
  
  if(!is.null(cb)) cb$clear()
  
  return(output)
}


# Summary method
summary.My_GD <- function(object) {
  return(object$trace)
}


# Print method
print.My_GD <- function(object){
  cat("Optimal parameters:\n")
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
    ggplot(GD_plot_df, aes(x = .time, y = value)) +
      geom_line() +
      geom_point(shape = 2) +
      scale_y_log10() +
      labs(title = "Negative marginal log likelihood vs Time", x = "Time", y = "Negative marginal log likelihood") +
      theme_bw()
    
  } else if (plot_no == 2) {
    
    par_norm_diff = sqrt((object$trace$par0.1 - object$trace$par1.1)^2 + (object$trace$par0.2 - object$trace$par1.2)^2)
    
    #browser()
    
    ggplot(object$trace, aes(x = i, y = par_norm_diff)) +
      geom_point() +
      scale_y_log10() +
      labs(title = "Norm of difference between estimated parameters vs Iteration", x = "Iteration", y = "Norm of difference") +
      theme_bw()
    
  } else if (plot_no == 3) {
    
    # Something with difference between last params and starting params
    par0 <- object$par0
    est <- object$est
    
    mu_grid <- seq(min(par0[1], est[1]) - 2, max(par0[1], est[1]) + 2, length.out = 50)
    sigma_grid <- seq(max(0.1, min(par0[2], est[2]) -2), max(par0[2], est[2]) + 2, length.out = 50)
    
    # Create a grid for mu and sigma
    grid <- expand.grid(mu = mu_grid, sigma = sigma_grid)
    
    #browser()
    # Calculate marginal_log_likelihood for each combination of mu and sigma
    grid$log_likelihood <- H_mult(x = x, par = cbind(grid$mu, grid$sigma), nu = nu)
    
    breaks_seq <- seq(min(grid$log_likelihood), max(grid$log_likelihood), length.out = 30)
    
    
    ggplot(grid) +
      geom_contour_filled(aes(x = mu, y = sigma, z = log_likelihood), breaks = breaks_seq, show.legend = FALSE) +
      geom_line(data = object$trace, aes(x = par0.1, y = par0.2), color = "darkblue", lwd = 1) +
      geom_point(data = object$trace, aes(x = par0.1, y = par0.2), color = "blue", shape = 1) +
      labs(title = "Contour Plot of Marginal Log-Likelihood", x = "Mu", y = "Sigma") +
      theme_bw()
    
  } else if ("true_par" %in% names(object$additional_args) & plot_no == 4) {
    
    true_par <- object$additional_args$true_par
    
    mu_oracle <- true_par[1]
    sigma_oracle <- true_par[2]
    
    par_norm_diff = sqrt((mu_oracle - object$trace$par0.1)^2 + (sigma_oracle - object$trace$par0.2)^2)
    
    ggplot(GD_plot_df, aes(x = .time, y = par_norm_diff)) +
      geom_line() +
      geom_point(shape = 2) +
      scale_y_log10() +
      labs(title = "Distance to true parameter vs Time", x = "Time", y = "Distance to true parameters") +
      theme_bw()
    
  } else {
    stop("Invalid plot number. Choose 1 2, 3 or 4 if you have supplied true_par")
  }
}
















