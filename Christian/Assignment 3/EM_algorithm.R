alpha_prime_func <- function(nu){
  return((nu + 1)/2)
}

beta_prime_func <- Vectorize(function(x, nu, mu, sigma){
  return(1/2 * (1 + (x - mu)^2 / (nu * sigma^2)))
}, vectorize.args = "x")



mu_hat_func <- function(x, beta_prime){
  return(sum(x / beta_prime) / sum(1 / beta_prime))
}

sigma_hat_func <- function(x, nu, mu_hat, alpha_prime, beta_prime){
  return(sqrt(alpha_prime/nu * mean((x - mu_hat)^2 / beta_prime)))
}




# EM_algorithms
EM_algorithm <- function(
    x, 
    par0, 
    nu, 
    maxiter = 1000, 
    tolerance = 1e-6, 
    cb = NULL){
  
  mu_hat <- par0[1]
  sigma_hat <- par0[2]
  
  alpha_prime <- alpha_prime_func(nu)
  
  for(i in 1:maxiter){
    mu_old <- mu_hat
    sigma_old <- sigma_hat
    beta_prime <- beta_prime_func(x = x, nu = nu, mu = mu_old, sigma = sigma_old)
    
    if(!is.null(cb)) cb$tracer()
    
    mu_hat <- mu_hat_func(x, beta_prime)
    sigma_hat <- sigma_hat_func(x, nu, mu_hat, alpha_prime, beta_prime)
    
    if(abs(mu_hat - mu_old) < tolerance & abs(sigma_hat - sigma_old) < tolerance){
      break
    }
    
    mu <- mu_hat
    sigma <- sigma_hat

  }
  
  return(c(mu, sigma))
}


EM_algorithm_cpp_naive <- function(
    x, 
    par0, 
    nu, 
    maxiter = 1000, 
    tolerance = 1e-6, 
    cb = NULL){
  
  mu_hat <- par0[1]
  sigma_hat <- par0[2]
  
  alpha_prime <- alpha_prime_func(nu)
  
  for(i in 1:maxiter){
    mu_old <- mu_hat
    sigma_old <- sigma_hat
    beta_prime <- beta_prime_func_cpp(x = x, nu = nu, mu = mu_old, sigma = sigma_old)
    
    if(!is.null(cb)) cb$tracer()
    
    mu_hat <- mu_hat_func(x, beta_prime)
    sigma_hat <- sigma_hat_func(x, nu, mu_hat, alpha_prime, beta_prime)
    
    if(abs(mu_hat - mu_old) < tolerance & abs(sigma_hat - sigma_old) < tolerance){
      break
    }
    
    mu <- mu_hat
    sigma <- sigma_hat
    
  }
  
  return(c(mu, sigma))
}



EM_algorithm_cpp <- function(
    x, 
    par0, 
    nu, 
    maxiter = 1000, 
    tolerance = 1e-6, 
    cb = NULL){
  
  mu_hat <- par0[1]
  sigma_hat <- par0[2]
  
  alpha_prime <- alpha_prime_func_cpp(nu)
  
  for(i in 1:maxiter){
    mu_old <- mu_hat
    sigma_old <- sigma_hat
    beta_prime <- beta_prime_func_cpp(x = x, nu = nu, mu = mu_old, sigma = sigma_old)
    
    
    mu_hat <- mu_hat_func_cpp(x, beta_prime)
    sigma_hat <- sigma_hat_func_cpp(x, nu, mu_hat, alpha_prime, beta_prime)
    
    if(!is.null(cb)) cb$tracer()
    
    if(abs(mu_hat - mu_old) < tolerance & abs(sigma_hat - sigma_old) < tolerance){
      break
    }
    
    mu <- mu_hat
    sigma <- sigma_hat
    
  }
  
  return(c(mu, sigma))
}





# Objects
EM_tracer <- tracer(c("mu_hat", "sigma_hat", "mu_old", "sigma_old", "i"), Delta = 0)

# EM class
EM <- function(x, 
               par0,
               nu, 
               maxiter = 1000, 
               tolerance = 1e-6, 
               cb = NULL, ...) {
  output <- structure(
    list(
      x = x,
      nu = nu,
      est = EM_algorithm_cpp(x = x, par0 = par0, nu = nu, maxiter = maxiter, tolerance = tolerance, cb = cb),
      trace = summary(cb),
      par0 = par0,
      additional_args = list(...)),
    class = "My_EM"
  )
  
  if(!is.null(cb)) cb$clear()
  
  return(output)
}


# Summary method
summary.My_EM <- function(object) {
  return(object$trace)
}


# Print method
print.My_EM <- function(object){
  cat("Optimal parameters:\n")
  print(object$est)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,3])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,4])
}




# Plot method
plot.My_EM <- function(object, plot_no = 1, ...) {
  x <- object$x
  nu <- object$nu
  
  
  neg_marginal_log_likelihood <- H_mult(x = x, par = object$trace[,1:2], nu = nu)
  
  
  EM_plot_df <- data.frame(object$trace, neg_marginal_log_likelihood)
  
  if (plot_no == 1) {
    ggplot(EM_plot_df, aes(x = .time, y = neg_marginal_log_likelihood)) +
      geom_line() +
      geom_point(shape = 2) +
      scale_y_log10() +
      labs(title = "Negative marginal log likelihood vs Time", x = "Time", y = "Negative marginal log likelihood") +
      theme_bw()
    
  } else if (plot_no == 2) {
    
    par_norm_diff = sqrt((object$trace$mu_old - object$trace$mu_hat)^2 + (object$trace$sigma_old - object$trace$sigma_hat)^2)
    
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
      geom_line(data = object$trace, aes(x = mu_old, y = sigma_old), color = "darkred", lwd = 1) +
      geom_point(data = object$trace, aes(x = mu_old, y = sigma_old), color = "red", shape = 1) +
      labs(title = "Contour Plot of Marginal Log-Likelihood", x = "Mu", y = "Sigma") +
      theme_bw()
    
  } else if ("true_par" %in% names(object$additional_args) & plot_no == 4) {
    
    true_par <- object$additional_args$true_par
    
    mu_oracle <- true_par[1]
    sigma_oracle <- true_par[2]
    
    par_norm_diff = sqrt((mu_oracle - object$trace$mu_hat)^2 + (sigma_oracle - object$trace$sigma_hat)^2)
    
    ggplot(EM_plot_df, aes(x = .time, y = par_norm_diff)) +
      geom_line() +
      geom_point(shape = 2) +
      scale_y_log10() +
      labs(title = "Distance to true parameter vs Time", x = "Time", y = "Distance to true parameters") +
      theme_bw()
    
  } else {
    stop("Invalid plot number. Choose 1 2, 3 or 4 if you have supplied true_par")
  }
}










