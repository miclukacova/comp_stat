## EM algorithm ##################################################

## E step

E.step <- function(X, par, nu){
  (nu + 1) / (1 + (X - par[1])^2 / (nu * par[2]))
}

## M step

M.step <- function(E_W, X, nu) {
  mu_k <- sum(E_W * X) / sum(E_W)
  sigma2_k <- sum(E_W * (X - mu_k)^2) / (length(X) * nu)
  return(c(mu_k, sigma2_k))
}


## EM algorithm

em_factory <- function(e_step, m_step, eps = 1e-6, nu) {
  force(e_step); force(m_step); force(eps); force(nu)
  function(par, X, epsilon = eps, cb = NULL, ...) {
    k <- 1
    repeat {
      if (!is.null(cb)) cb()
      k <- k + 1
      par0 <- par
      E_W <- e_step(X = X, par = par, nu = nu)
      par <- m_step(E_W = E_W, X = X, nu = nu)
      if (sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon))
        break
    }
    par  # Returns the parameter estimate
  }
}


# Tracer
EM_tracer <- tracer(c("par", "par0"), Delta = 0) 



## S3 class

EM <- function(par0, X, eps = 1e-6, 
               M_step = M.step,
               E_step = E.step,
               par_true,
               cb = EM_tracer, 
               ...) {
  
  EM_alg <- em_factory(e_step = E_step, m_step = M_step, eps = 1e-6, nu = par_true[3])
  par_est <- EM_alg(par = par0, X = X, cb = cb$tracer, ...)
  
  EM_trace <- summary(cb)
  
  EM_trace <- transform(
    EM_trace,
    n = seq_len(nrow(EM_trace)),
    par_norm_diff = sqrt((par0.1 - par.1)^2 + (par0.2 - par.2)^2)
  )
  
  cb$clear()
  
  structure(
    list(
      est = par_est,
      algorithm = EM_alg,
      trace = EM_trace,
      start_par = par0,
      par_true = par_true,
      X = X,
      additional_args = list(...)),
    class = "My_EM"
  )
}

print.My_EM <- function(x) {
  cat("EM algorithm: \n")
  cat("Estimated parameters: \n")
  cat("mu: ", x$est[1], "\n")
  cat("sigma2: ", x$est[2], "\n")
  cat("True parameters: \n")
  cat("mu: ", x$par_true[1], "\n")
  cat("sigma2: ", x$par_true[2], "\n")
  cat("nu: ", x$par_true[3], "\n")
}





# Plot method
plot.My_EM <- function(object, plot_no = 1, mle_par = NULL,...) {
  
  plot_data <- object$trace %>%
    mutate("Dist_mu" = abs(par.1 - object$par_true[1]),
           "Dist_sigma2" = abs(par.2 - object$par_true[2])) 
  
  plot_data$loglikelihood <- apply(cbind(plot_data$par.1, plot_data$par.2), 1, function(params) {
    loglik(x = object$X, params, nu = object$par_true[3])})
  
  plot_data$dist_loglik <- abs(plot_data$loglikelihood - loglik(x = object$X, mle_par, nu = object$par_true[3]))
  
  if(!is.null(mle_par)) {
    plot_data <- plot_data %>%
      mutate("Dist_mle_mu" = abs(par.1 - mle_par[1]),
             "Dist_mle_sigma2" = abs(par.2 - mle_par[2]))
    
  }

  if (plot_no == 1) {
    p <- ggplot(plot_data, aes(x = .time)) +
      geom_line(aes(y = Dist_mu, color = "mu")) +
      geom_line(aes(y = Dist_sigma2, color = "sigma2")) +
      scale_color_manual(values = c("mu" = "red3", "sigma2" = "blue3")) +
      scale_y_log10() +
      labs(title = "Dist to true par vs Time", x = "Time", y = "Dist")
  }
  
  if (plot_no == 2) {
    p <- ggplot(plot_data, aes(x = .time)) +
      geom_line(aes(y = Dist_mle_mu, color = "mu")) +
      geom_line(aes(y = Dist_mle_sigma2, color = "sigma2")) +
      scale_color_manual(values = c("mu" = "red3", "sigma2" = "blue3")) +
      scale_y_log10() +
      labs(title = "Dist to MLE vs Time", x = "Time", y = "Dist")
  }
  
  if (plot_no == 3) {
    p <- ggplot(plot_data, aes(x = .time)) +
      geom_line(aes(y = dist_loglik), color = "blue4") +
      scale_y_log10() +
      labs(title = "Loglik - Loglik MLE vs. Time", x = "Time", y = "Abs. Dist")
  }
  
  if (plot_no == 4) {
    p <- ggplot(plot_data, aes(x = .time)) +
      geom_line(aes(y = loglikelihood), color = "blue4") +
      labs(title = "Log likelihood vs. Time", x = "Time", y = "Log likelihood")
  }
  if (plot_no == 5) {
    
    plot_data$neg_loglik <- -plot_data$loglikelihood * 1/length(object$X)
    
    p <- ggplot(plot_data, aes(x = .time)) +
      geom_line(aes(y = neg_loglik), color = "blue4") +
      scale_y_log10() +
      labs(title = "Neg. Loglikelihood vs. Time", x = "Time", y = "Neg. Loglik")
  }
  
  return(p)
}


# Method to extract plot data
plot_data <- function(x) {
  UseMethod("plot_data")
}

plot_data.My_EM <- function(object, mle_par = NULL) {
  
  plot_data <- object$trace %>%
    mutate("Dist_mu" = abs(par.1 - object$par_true[1]),
           "Dist_sigma2" = abs(par.2 - object$par_true[2])) 
  
  plot_data$loglikelihood <- apply(cbind(plot_data$par.1, plot_data$par.2), 1, function(params) {
    loglik(x = object$X, params, nu = object$par_true[3])})
  
  plot_data$neg_loglik <- -plot_data$loglikelihood * 1/length(object$X)
  
  if(!is.null(mle_par)) {
    plot_data <- plot_data %>%
      mutate("Dist_mle_mu" = abs(par.1 - mle_par[1]),
             "Dist_mle_sigma2" = abs(par.2 - mle_par[2]))
    
    plot_data$dist_loglik <- abs(plot_data$loglikelihood - loglik(x = object$X, mle_par, nu = object$par_true[3]))
  }
  
  return(plot_data)
}


# Heatmap method

heat_map <- function(obj, mle, log_lik = loglik, x, grid_vals_m = c(0, 3), grid_vals_s = c(1, 4)) {
  UseMethod("heat_map")
}

heat_map.default <- function(obj, ...) {
  stop("heat_map is not implemented for objects of class ", class(obj))
}

heat_map.My_EM <- function(obj, mle, log_lik = loglik, x, grid_vals_m = c(0, 3), grid_vals_s = c(1, 4)){
  
  # Create a grid of m and s values
  m_values <- seq(grid_vals_m[1], grid_vals_m[2], length.out = 100)
  s_values <- seq(grid_vals_s[1], grid_vals_s[2], length.out = 100)
  
  # Create a dataframe to store the values of m, s, and loglik
  results <- expand.grid(m = m_values, s = s_values)
  results$log_lik <- apply(results, 1, function(row) log_lik(x = x, row, 1))
  
  # Path of EM
  em_path <- obj$trace %>% select(c(par.1,par.2)) %>% rename(mu = par.1, sigma = par.2)
  em_path$log_lik <- apply(em_path, 1, function(row) log_lik(x = x, row, 1))
  

  
  # Plot the heatmap with contours and a point using ggplot2
  ggplot(results, aes(x = m, y = s, fill = log_lik)) +
    geom_tile() +
    geom_contour(aes(z = log_lik), color = "darkseagreen4", bins = 50) +
    scale_fill_gradient2(mid = "darkseagreen1", high = "coral1", low = "seagreen4", midpoint = -3500) +
    geom_point(aes(x = mle[1], y = mle[2], colour = "Full data MLE"), size = 2.5) +
    geom_point(aes(x = obj$est[1], y = obj$est[2], colour = "EM"), size = 2.5) +
    geom_path(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
    geom_point(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
    scale_color_manual(values = c("EM" = "seagreen", "Full data MLE" = "seagreen3")) +
    labs(x = "mu", y = "sigma", fill = "loglikelihood") +
    theme_minimal()
}

