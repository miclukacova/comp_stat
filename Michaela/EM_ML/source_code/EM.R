## EM algorithm ##################################################

## E step

E.step <- function(X, par, nu){
  (nu + 1) / (1 + (X - par[1])^2 / (nu * par[2]))
}

## M step

M.step <- function(E.W, X, nu) {
  mu_k <- sum(E.W * X) / sum(E.W)
  sigma2_k <- sum(E.W * (X - mu_k)^2) / (length(X) * nu)
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
      E.W <- e_step(X = X, par = par, nu = nu)
      par <- m_step(E.W = E.W, X = X, nu = nu)
      if (!is.null(cb)) cb()
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
    par_norm_diff = sqrt((par.1 - par0.1)^2 + (par.2 - par0.2)^2)
  )
  
  output = structure(
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
  cb$clear()
  return(output)
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
