test_estimator <- function(estimator, true_params = NULL) {
  
  if(is.null(true_params)){
    true_params <- data.frame(mu = c(1, 5, 0), sigma2 = c(3, 10, 3), nu = c(1,1,3))
  }
  
  # Create empty vectors to store estimated values
  mu_est <- numeric(6*3) # rows are parameter index, columns are sample size
  sigma2_est <- numeric(6*3) # rows are parameter index, columns are sample size
  
  sample_size = rep(c(20, 30, 40, 50, 100, 200), each = 3)
  param_id = rep(1:3, times = 6)
  
  param_labels <- c(
    "1" = "(1,3,1)",
    "2" = "(5,10,1)",
    "3" = "(0,3,3)"
  )
  
  for(i in 1:(length(sample_size)) ? nrow(true_params)){
    ss <- sample_size[i]
    par <- true_params[param_id[i],]
    data <- simulate(n = ss, par = par) 
    mles <- my_mle(data, nu = par$nu)
    mu_est[i] <- mles$mu_mle
    sigma2_est[i] <- mles$sigma2_mle
  }
  
  data <- data.frame(
    sample_size = sample_size,
    param_id = param_id,  # 3 parameters, each with 6 sample sizes
    mu_mle = mu_est,
    sigma2_mle = sigma2_est,
    mu_true = true_params[param_id,1],  # True mu value
    sigma2_true = true_params[param_id,2]  # True sigma value
  )
  
  # Reshape data to long format
  data_long <- data %>%
    pivot_longer(cols = starts_with("mu_") | starts_with("sigma2_"), 
                 names_to = c("Parameter", "Type"), 
                 names_sep = "_", 
                 values_to = "Value")
  
  ggplot(data_long, aes(x = sample_size, y = Value, color = Type)) +
    geom_line(aes(linetype = Type), size = 1) +
    geom_point(size = 2) +
    facet_grid(Parameter ~ param_id, scales = "free_y", 
               labeller = labeller(param_id = param_labels)) +  # Separate rows for `mu` and `sigma`
    scale_color_manual(values = c("mle" = "blue", "true" = "red")) +
    scale_linetype_manual(values = c("mle" = "solid", "true" = "dashed")) +
    labs(x = "Sample Size", y = "Parameter Value", 
         title = "MLE vs. True Values for Parameters Across Sample Sizes") +
    theme_minimal() +
    theme(legend.position = "top", legend.title = element_blank())
  
}