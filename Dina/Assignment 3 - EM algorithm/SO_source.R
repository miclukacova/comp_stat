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

## EM ###############################################

EM_alg <- function(x, param, ny, max_iter = 20, epsilon = 1e-10, cb = NULL){
  mu_mark <- param$mu #Getting initial parameters
  sigma_mark <- param$sigma
  ny <- ny
  k <- (ny + 1) / 2   #Defining i and k, which only depends on ny:
  i <- 0
  
  while (i < max_iter) {
    mu_old <- mu_mark
    sigma_old <- sigma_mark
    t_old <- 2 / (1 + (x - mu_old)^2 / (ny * sigma_old))
    
    mu_mark <- sum(t_old * x) / sum(t_old)
    sigma_mark <- 1 / (length(x) * ny) * k * sum(t_old * (x - mu_mark)^2)
    #Calling  cb
    if(!is.null(cb)) cb()
    #Checking convergence
    if (sum((c(mu_mark, sigma_mark) - c(mu_old, sigma_old))^2) 
        <= epsilon * (sum(c(mu_mark, sigma_mark)^2) + epsilon)) break
    i <- i + 1
  }
  c(mu_mark, sigma_mark)
}

## log-likelihood  ###########################################

log_lik <- function(m, s, ny, x = samp$x){
  k <- (ny + 1) / 2
  return( - length(x) * log(sqrt(s)) - sum(k * log(1 + (x - m)^2 / (ny * s))))
}

## Norm of the difference between parameters ###########################################

par_norm_diff <- function(mu_old, sigma_old, mu_mark, sigma_mark){
  return(sqrt((mu_old - mu_mark)^2 + (sigma_old - sigma_mark)^2))
}

## Gradient  ###########################################

grad_old <- function(m, s, ny, x){
  k <- (ny + 1) / 2
  t <- mean(2 / (1 + (x - m)^2 / (ny * s)))
  grad_m <- - t * mean(k * (x - m) / (ny * s))
  grad_s <- - mean(1 / s) + t * mean(k * (x - m)^2 / (ny * sqrt(s)^3))
  return(c(grad_m, grad_s))
}

grad <- function(m, s, ny, x = samp$x) {
  n <- length(x)
  k <- (ny + 1) / 2
  
  # Partial derivative with respect to m
  partial_m <- -sum((k * 2 * (x - m)) / (ny * s * (1 + (x - m)^2 / (ny * s))))
  
  # Partial derivative with respect to s
  partial_s <- -n / (2 * s) - sum(k * (x - m)^2 / (ny * s^2 * (1 + (x - m)^2 / (ny * s))))
  
  return(c(partial_m, partial_s))
}


## Sampler  ##############################################

samp_our_distribution <- function(n, param){
    w <- rchisq(n, df = param$ny)
    x <- rnorm(n, mean = param$mu, sd = sqrt((param$ny * param$sigma) / w))
    y <- cbind(x, w)
    return(data.frame(y))
  }

## "Theoretical parameters" #####################################

theo_par <- function(x, w, ny){
  mu_opt <- sum(w * x) / sum(w)
  sigma_opt <- 1/(length(x) * ny) * sum(w * (x - mu_opt)^2)
  return(data.frame(mu = mu_opt, sigma = sigma_opt ))
}

##### Fischer Information #####################################




##### Tracer #####################################

EM_tracer <- tracer(c("mu_mark", "sigma_mark", "mu_old", "sigma_old", "ny", "i"), Delta = 0)

log_lik_mult <- function(m, s, ny, x){
  param <- cbind(m, s, ny)
  print(param)
  apply(param, 1, function(par) log_lik(par[1], par[2], par[3], x = x))
}

grad_mult <- function(m, s, ny, x){
  param <- cbind(m, s, ny)
  apply(param, 1, function(par) grad(par[1], par[2], par[3], x = x))
}

### Suboptimiality ###########################################


###### S3 Classes ###########################################

# Parameters class

parameters <- function(mu, sigma, ny) {
  structure(
    list(
      mu = mu,
      sigma = sigma,
      ny = ny, 
      par = c(mu, sigma, ny)),
    class = "parameters"
  )
}

sim <- function(x) {
  UseMethod("sim")
}

sim <- function(object, N) {
  return(as.data.frame(samp_our_distribution(n = N, param = object)))
}

sim(parameters(1, 2, 3), 100)
#-------------------------------------------------------------------------------
#                                   EM class                                  #
#-------------------------------------------------------------------------------

# EM class
EM <- function(x, par0, max_iter, epsilon, cb = NULL, par_true,...) {
  structure(
    list(
      est = EM_alg(x = x, param = par0, max_iter = max_iter, epsilon = epsilon, cb = cb$tracer,...),
      trace = summary(cb),
      start_par = par0,
      x = x,
      additional_args = list(...)),
    class = "My_EM"
  )
}


# EM class
EM <- function(x, par0, max_iter, epsilon, cb = EM_tracer, par_true,...) {
  
  
  par_est <- EM_alg(x = x, param = par0, ny = par_true$ny, max_iter = max_iter, epsilon = epsilon, cb = cb$tracer)
  
  EM_trace <- summary(cb)
  
  EM_trace <- transform(
    EM_trace,
    par_norm_diff = sqrt((mu_old - mu_mark)^2 + (sigma_old - sigma_mark)^2),
    loglik = log_lik(mu_mark, sigma_mark, ny = par_true$ny, x = x)
  )
  
  output = structure(
    list(
      est = par_est,
      trace = EM_trace,
      start_par = par0,
      par_true = par_true,
      x = x,
      additional_args = list(...)),
    class = "My_EM"
  )
  cb$clear()
  return(output)
}


#test <- EM(sim1$x, parameters(1, 2, 3), 20, 1e-10, par_true = parameters(1, 1,1))

#test$trace
# Summary method
summary.My_EM <- function(object) {
  return(object$trace[c(1,2,5,6,7)])
}


# Print method
print.My_EM <- function(object){
  cat("EM algorithm:\n")
  cat("True parameters:\n")
  print(object$par_true$par)
  cat("\nOptimal parameters:\n")
  print(object$est)
  cat("\nNumber of iterations:\n")
  print(tail(object$trace, 1)[,6])
  cat("\nTotal time:\n")
  print(tail(object$trace, 1)[,7])
}

base <- "#4f7942"
my_theme <-   theme(
  text = element_text(size = 16),  # Change the base text size
  plot.title = element_text(size = 18),  # Title size
  axis.title = element_text(size = 16),  # Axis titles size
  axis.text = element_text(size = 14),  # Axis text size
  legend.title = element_text(size = 18),  # Legend title size
  legend.text = element_text(size = 16),  # Legend text size
  line = element_line(linewidth = 4) +
    theme_bw()
)
#Plot method for convergence plots
plot.My_EM <- function(object, plot_no = 1, mle_est = NULL,...) {

  if(!is.null(mle_est)) {
    plot_data <- object$trace %>%
      mutate("Dist_mle_mu" = abs(mu_mark - mle_est$mu),
             "Dist_mle_sigma2" = abs(mu_mark - mle_est$sigma),
             "suboptimality" = abs(loglik - log_lik(mle_est$mu, mle_est$sigma, ny = object$par_true$ny, x = object$x)))
  }

  if (plot_no == 1) {
    p <- ggplot(data = object$trace, aes(x = .time, y = par_norm_diff)) +
      geom_line(color = "coral4", size = 1.2) +
      geom_point(color = "coral4", size = 2.5) +
      labs(
        x = "Time",
        y = expression(log(bgroup("||", theta^{(n-1)} - theta^n, "||")))
      ) +
      scale_y_log10() +
      my_theme
    return(p)
  }
  if (plot_no == 4) {
    p <- ggplot(data = object$trace, aes(x = i, y = par_norm_diff)) +
      geom_line(color = base, size = 1.2) +
      geom_point(color = base, size = 2.5) +
      labs(
        x = "Iterations",
        y = expression(log(bgroup("||", theta^{(n-1)} - theta^n, "||")))
      )+
      scale_y_log10() +
      my_theme
    return(p)
  }
  else if (plot_no == 2) {
    p <- ggplot(data = object$trace, aes(x = .time, y = -loglik)) +
      geom_line(color = base, size = 1.2) +
      geom_point(color = base, size = 2.5) +
      labs(
           x = "Time",
           y = "Log-likelihood - logscale") +
      scale_y_log10() +
      my_theme
    return(p)
  }
  else if (plot_no == 3) {
    p <- ggplot(data = plot_data, aes(x = .time, y = suboptimality)) +
      geom_line(color = "steelblue4", size = 1.2) +
      geom_point(color = "steelblue4", size = 2) +
      labs(
           x = "Time",
           y = "Suboptimality - logscale") +
      scale_y_log10() +
      my_theme
    return(p)
  }
  else {
    stop("Invalid plot number")
  }
}
#plot(test, 3, mle_est = theo_par(test$x, sim1$w, test$par_true$ny))

heatmap <- function(x) {
  UseMethod("heatmap")
}

grad <- function(mu, sigma, nu, x){
  n <- length(x)
  d_mu <- - (nu + 1) * sum((x - mu) / (nu * sigma^2 + (x - mu)^2))
  d_sigma <- n/sigma - (nu + 1) / sigma * sum((x - mu)^2 / (nu * sigma^2 + (x - mu)^2))
  return(c(d_mu, d_sigma))
}

heatmap.My_EM <- function(object, mle_est, gradient = F, path = F, ...){
  
  # Create a grid of m and s values
  m_values <- seq(-3, 3, length.out = 100)
  s_values <- seq(0.01, 5.5, length.out = 100)
  
  # Create a dataframe to store the values of m, s, and log_lik
  results <- expand.grid(m = m_values, s = s_values)
  results$log_lik <- apply(results, 1, function(row) log_lik(row['m'], row['s'],  ny = object$par_true$ny, x = as.matrix(object$x)))
  EM_param <- object$est
  theo_param <- mle_est 
  
  if(gradient == F & path ==F){
  # Plot the heatmap with contours and a point using ggplot2
  p <- ggplot(results, aes(x = m, y = s, fill = log_lik)) +
    geom_tile() +
    geom_contour(aes(z = log_lik), color = "darkseagreen4", bins = 50) +
    scale_fill_gradient2(mid = "darkseagreen1", high = "coral1", low = "seagreen4", midpoint = -3000) +
    geom_point(aes(x = theo_param$mu, y = theo_param$sigma), colour = "seagreen", size = 2.5) +
    geom_point(aes(x = EM_param[1], y = EM_param[2]), colour = "coral4", size = 2.5) +
    # Add text labels next to the points
    geom_text(aes(x = theo_param$mu, y = theo_param$sigma, label = "True MLE"), 
              hjust = 1.2, vjust = 0.5, colour = "seagreen", size = 4) +
    geom_text(aes(x = EM_param[1], y = EM_param[2], label = "EM Estimate"), 
              hjust = -0.2, vjust = 0.5, colour = "coral4", size = 4) +
    labs(x = expression(mu), y = expression(sigma), fill = "loglikelihood") +
    theme_minimal()
  
  return(p)
  }
  else if(gradient == T & path == F){
    
    
    # Create a dense grid for the heatmap (100 x 100)
    m_values <- seq(-3, 3, length.out = 100)
    s_values <- seq(0.01, 5.5, length.out = 100)
    results <- expand.grid(m = m_values, s = s_values)
    results$log_lik <- apply(results, 1, function(row) log_lik(row['m'], row['s'], ny = object$par_true$ny, x = as.matrix(object$x)))
    
    # Create a coarser grid for gradient arrows (10 x 10)
    m_values_coarse <- seq(-3, 3, length.out = 10)
    s_values_coarse <- seq(0.01, 5.5, length.out = 10)
    coarse_grid <- expand.grid(m = m_values_coarse, s = s_values_coarse)
    
    # Calculate gradients on the coarse grid
    gradients <- grad_mult(coarse_grid$m, coarse_grid$s, ny = object$par_true$ny, x = as.matrix(object$x))
    
    # Calculate the lengths of the gradient vectors
    gradient_lengths <- -sqrt(gradients[1, ]^2 + gradients[2, ]^2)
    
    # Replace zeros in gradient_lengths with a small value to avoid division by zero
    gradient_lengths[gradient_lengths == 0] <- 1e-6  # Small value
    
    # Normalize gradients to have a fixed length of approximately 0.1
    arrow_length <- 0.3  # Desired length of arrows
    normalized_gradients <- gradients / matrix(gradient_lengths, nrow = 2, ncol = length(gradient_lengths), byrow = TRUE) * arrow_length
    
    # Calculate the end points of the arrows based on the normalized gradients
    coarse_grid <- transform(coarse_grid,
                             m_end = m + normalized_gradients[1, ],  # Add normalized gradient in the m direction
                             s_end = s + normalized_gradients[2, ])  # Add normalized gradient in the s direction
    
    # Get EM and true parameter estimates
    EM_param <- object$est
    theo_param <- mle_est
    
    # Plot the heatmap with contours, points, and gradient arrows using ggplot2
    p <- ggplot(results, aes(x = m, y = s)) +
      geom_tile(aes(fill = log_lik)) +  # Specify `fill` only here
      geom_contour(aes(z = log_lik), color = "darkseagreen4", bins = 50) +
      scale_fill_gradient2(mid = "darkseagreen1", high = "coral1", low = "seagreen4", midpoint = -3000) +
      geom_point(aes(x = theo_param$mu, y = theo_param$sigma), colour = "darkgreen", size = 5) +
      geom_point(aes(x = EM_param[1], y = EM_param[2]), colour = "coral4", size = 5) +
      
      # Add text labels next to the points
      geom_text(aes(x = theo_param$mu, y = theo_param$sigma, label = "True MLE"), 
                hjust = 1.2, vjust = 0.5, colour = "darkgreen", size = 5) +
      geom_text(aes(x = EM_param[1], y = EM_param[2], label = "EM Estimate"), 
                hjust = -0.2, vjust = 0.5, colour = "coral4", size = 5) +
      
      # Add arrows representing the gradient on the coarse grid
      geom_segment(data = coarse_grid,
                   aes(x = m, y = s, xend = m_end, yend = s_end),
                   arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.7) +
      
      labs(x = expression(mu), y = expression(sigma), fill = "loglikelihood") +
      theme_minimal()
    
    return(p)
  }
  else if(gradient == F & path == T){
    # Path of EM
    em_path <- object$trace %>% select(c(mu_old,sigma_old)) %>% rename(mu = mu_old, sigma = sigma_old)
    em_path$log_lik <- log_lik(m = em_path$mu, s = em_path$sigma,  ny = object$par_true$ny, x = as.matrix(object$x))
    
    p <- ggplot(results, aes(x = m, y = s, fill = log_lik)) +
      geom_tile() +
      geom_contour(aes(z = log_lik), color = "darkseagreen4", bins = 50) +
      scale_fill_gradient2(mid = "darkseagreen1", high = "coral1", low = "seagreen4", midpoint = -2700) +
      #geom_path(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
      geom_point(data = em_path, aes(x = mu, y = sigma),color = "seagreen", size = 3, ) +
      geom_point(aes(x = theo_param$mu, y = theo_param$sigma),colour = "darkgreen", size = 5, alpha = 0.6) +
      geom_point(aes(x = EM_param[1], y = EM_param[2]), colour = "coral4", size = 5) +
      # Add text labels next to the points
      geom_text(aes(x = theo_param$mu, y = theo_param$sigma, label = "True MLE"), 
                hjust = 1.2, vjust = 0.5, colour = "darkgreen", size = 5, alpha = 0.6) +
      geom_text(aes(x = EM_param[1], y = EM_param[2], label = "EM Estimate"), 
                hjust = -0.2, vjust = 0.5, colour = "coral4", size = 5) +
      labs(x = expression(mu), y = expression(sigma), fill = "loglikelihood") +
      xlim(c(-3,3)) +
      ylim(c(0,5))
      theme_minimal()
    return(p)

  }
}


#test$trace
#heatmap.My_EM(test, theo_par(test$x, sim1$w, test$par_true$ny), gradient = T)
#grad_mult(theo_par(test$x, sim1$w, 1)$mu, theo_par(test$x, sim1$w, 1)$sigma, 3, test$x)
# Method to extract plot data
plot_data <- function(x) {
  UseMethod("plot_data")
}


plot_data.My_SGD <- function(object) {
  x <- object$additional_args$x
  y <- object$additional_args$y
  
  loss <- squared_error_mult(x = x, y = y, 
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


