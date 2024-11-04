### Kontur plots 
# Snak med Dina om dem her


H_contur <- function(alpha, beta, x, y) {
  mean((y - f_contur(alpha, beta, x))^2)
} 

f_contur <- function(alpha, beta, x) {
  gamma <- 1
  rho <- 2
  gamma + (rho - gamma)/(1 + exp(beta * log(x) - alpha))
}

# Create a grid of m and s values
a_values <- seq(-10, 10, length.out = 130)
b_values <- seq(-1, 2, length.out = 130)

# Create a dataframe to store the values of m, s, and loglik
results <- expand.grid(a = a_values, b = b_values)
results$loss <- apply(results, 1, function(row) H_contur(x = x, y = y, alpha = row[2], beta = row[2]))

# Path of EM
#em_path <- obj$trace %>% select(c(par.1,par.2)) %>% rename(mu = par.1, sigma = par.2)
#em_path$log_lik <- apply(em_path, 1, function(row) log_lik(x = x, row, 1))


# Plot the heatmap with contours and a point using ggplot2
p1 <- ggplot(results, aes(x = a, y = b, fill = loss)) +
  geom_tile() +
  geom_contour(aes(z = loss), color = "#10B9F1", bins = 20) +
  scale_fill_gradient2( low = "#132B43", high = "#70B1F7", space = "Lab",
                        na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  #geom_point(aes(x = mle[1], y = mle[2], colour = "Full data MLE"), size = 2.5) +
  #geom_point(aes(x = obj$est[1], y = obj$est[2], colour = "EM"), size = 2.5) +
  #geom_path(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
  #geom_point(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
  #scale_color_manual(values = c("EM" = "seagreen", "Full data MLE" = "seagreen3")) +
  labs(x = "alpha", y = "beta", fill = "Loss") +
  theme_minimal()


H_contur1 <- function(gamma, rho, x, y) {
  mean((y - f_contur1(gamma, rho, x))^2)
} 

f_contur1 <- function(gamma, rho, x) {
  alpha <- 2
  beta <- 5
  gamma + (rho - gamma)/(1 + exp(beta * log(x) - alpha))
}

# Create a grid of m and s values
g_values <- seq(0, 2, length.out = 130)
r_values <- seq(0, 4, length.out = 130)

# Create a dataframe to store the values of m, s, and loglik
results <- expand.grid(g = g_values, r = r_values)
results$loss <- apply(results, 1, function(row) H_contur1(x = x, y = y, gamma = row[2], rho = row[2]))

# Path of EM
#em_path <- obj$trace %>% select(c(par.1,par.2)) %>% rename(mu = par.1, sigma = par.2)
#em_path$log_lik <- apply(em_path, 1, function(row) log_lik(x = x, row, 1))


# Plot the heatmap with contours and a point using ggplot2
p2 <- ggplot(results, aes(x = g, y = r, fill = loss)) +
  geom_tile() +
  geom_contour(aes(z = loss), color = "#10B9F1", bins = 20) +
  scale_fill_gradient2( low = "#132B43", high = "#70B1F7", space = "Lab",
                        na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  #geom_point(aes(x = mle[1], y = mle[2], colour = "Full data MLE"), size = 2.5) +
  #geom_point(aes(x = obj$est[1], y = obj$est[2], colour = "EM"), size = 2.5) +
  #geom_path(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
  #geom_point(data = em_path, aes(x = mu, y = sigma), color = "seagreen", size = 1) +
  #scale_color_manual(values = c("EM" = "seagreen", "Full data MLE" = "seagreen3")) +
  labs(x = "gamma", y = "rho", fill = "Loss") +
  theme_minimal()