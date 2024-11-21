#Loglikelihood 
loglik <- function(x, par, nu) {
  mu <- par[1]
  sigma2 <- par[2]
  sum(- log(sqrt(sigma2)) - (nu + 1) / 2 * log(1 + (x - mu)^2 / (nu * sigma2)))
}

#loglik(sim1$x, c(1, 2), 1)
# Gradient of negative log likelihood
grad_negloglik <- function(mu, sigma, nu, x){
  n <- length(x)
  d_mu <- - (nu + 1) * sum((x - mu) / (nu * sigma^2 + (x - mu)^2))
  d_sigma <- n/sigma - (nu + 1) / sigma * sum((x - mu)^2 / (nu * sigma^2 + (x - mu)^2))
  return(c(d_mu, d_sigma))
}

#grad_negloglik(1, 2, 1, sim1$x)
# Negative loglikelihood
neg_loglik <- function(par, x, nu) {
  - loglik(par = par, x = x, nu)
}

#neg_loglik(c(1, 2), sim1$x, 1)
# Gradient descent algorithm 
grad_desc_mom <- function(
    par,
    grad = grad_negloglik,
    H = neg_loglik,
    t0 = 1,
    maxit = 1200,
    cb = NULL,
    epsilon = 1e-6,
    beta = 0.8,
    alpha = 0.1,
    x,
    nu,
    mu,
    ...) {
  
  n <- length(x)
  par_old <- par
  
  for (i in 1:maxit) {
    
    #browser()
    
    # Calculations of objective and gradient
    value <- 1 / n * H(x = x, par = par, nu = nu)
    gr <- 1/n * grad(mu = par[1], sigma = sqrt(par[2]), nu = nu, x = x)
    
    grad_norm <- sum(gr^2)
    
    t <- t0
    # Proposed descent step
    par_new <- par - t * gr + mu * (par - par_old)
    
    # Backtracking line search
    while (1 / n * H(x = x, par = par_new, nu = nu) > value - alpha * t * grad_norm) {
      t <- beta * t
      par_new <- par - t * gr
    }
    
    # Callback
    if (!is.null(cb)) cb()
    
    # Convergence criterion 
    if (sum((par_new - par)^2) <= epsilon * (sum(par_new^2) + epsilon)) break
    
    par_old <- par
    par <- par_new
  }
  
  if (i == maxit)  warning("Maximal number, ", maxit, ", of iterations reached")
  
  par
}

#grad_desc(par = c(1, 2), x = sim1$x, nu = 1)
##### Tracer #####################################

GD_tracer <- tracer(c("par_new", "par", "value", "gr", "grad_norm", "i"), Delta = 0) 


###### GD class ###########################################

GD_m <- function(par,
               grad = grad_negloglik,
               H = neg_loglik,
               t0 = 1,
               maxit = 1200,
               cb = GD_tracer,
               epsilon = 1e-6,
               beta = 0.8,
               alpha = 0.1,
               nu = NULL,
               mu = 0.9,
               ...) {
  
  est <- grad_desc_mom(par = par, 
                   grad = grad, 
                   t0 = t0,
                   maxit = maxit,
                   cb = cb$tracer,
                   epsilon = epsilon,
                   H = H,
                   beta = beta,
                   alpha = alpha,
                   nu = nu,
                   mu = mu,
                   ...)
  
  GD_trace = summary(GD_tracer)
  
  GD_trace <- transform(
    GD_trace,
    n = seq_len(nrow(GD_trace)),
    par_norm_diff = sqrt((par.1 - par_new.1)^2 + (par.2 - par_new.2)^2)
  )
  
  GD_tracer$clear()
  
  structure(
    list(
      est = est,
      trace = GD_trace,
      start_par = par,
      nu = nu,
      additional_args = list(...)),
    class = "My_GD_m"
  )
}


# Print method
print.My_GD_m <- function(object){
  cat("Estimated parameters:\n")
  print(object$est)
  cat("Number of iterations:\n")
  print(tail(object$trace, 1)[,9])
  cat("Total time:\n")
  print(tail(object$trace, 1)[,10])
}

base <- "#4f7942"
# Plot method
plot.My_GD_m <- function(object, plot_no = 1, mle_est = NULL,...) {
  x <- object$x
  nu <- object$nu
  GD_plot_df <- data.frame(object$trace)
  
  if(!is.null(mle_est)){
    GD_plot_df$abs_dist_from_par <- apply(GD_plot_df[,1:2], 1, FUN = function(par_est) sum(abs(par_est - mle_est)))
  }
  
  
  if (plot_no == 1) {
    ggplot() +
      geom_line(data = GD_plot_df, aes(x = .time, y = value), color = base, size = 1.2) +
      geom_point(data = GD_plot_df, aes(x = .time, y = value), color = base, size = 2.5) +
      scale_y_log10() +
      labs(
        x = "Time",
        y = "Log-likelihood - logscale")
  } else if (plot_no == 2){
    ggplot() +
      geom_line(data = GD_plot_df, aes(x = .time, y = par_norm_diff), color = "coral4", size = 1.2) +
      geom_point(data = GD_plot_df, aes(x = .time, y = par_norm_diff), color = "coral4", size = 2.5) +
      scale_y_log10() +
      labs(
        x = "Time",
        y = expression(log(bgroup("||", theta^{(n-1)} - theta^n, "||")))
      )
  } else if (plot_no == 3){
    ggplot() +
      geom_line(data = GD_plot_df, aes(x = .time, y = abs_dist_from_par), color = "steelblue4", size = 1.2) +
      geom_point(data = GD_plot_df, aes(x = .time, y = abs_dist_from_par), color = "steelblue4", size = 2.5) +
      scale_y_log10() +
      labs(
        x = "Time",
        y = "Suboptimality - logscale")
  } else {
    stop("Invalid plot number. Please choose 1 2 or 3.")
  }
}

# Method to extract plot data
plot_data <- function(x) {
  UseMethod("plot_data")
}


plot_data.My_GD_m <- function(object) {
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


