

source("~/Desktop/Studie/Master/First year/Block 1/Computational Statistics/comp_stat/Christian/Assignment 3/Generate_data.R", echo=TRUE)
source("~/Desktop/Studie/Master/First year/Block 1/Computational Statistics/comp_stat/Christian/Assignment 3/EM_algorithm.R", echo=TRUE)
source("~/Desktop/Studie/Master/First year/Block 1/Computational Statistics/comp_stat/Christian/Assignment 3/Gradient_descent.R", echo=TRUE)


################################################################################
#                           Generate_data.R tests                              #
################################################################################

N <- 1000000
nu <- 3
mean_oracle <- 7
sigma_oracle <- 5


oracle_param <- parameters(mu = mean_oracle, sigma = sigma_oracle, nu = nu)
oracle_par <- oracle_param$par
x <- simulate(oracle_param, N = N)$x
w <- simulate(oracle_param, N = N)$w

mean(x)
sd(x)
hist(x)

mean(w)
sd(w)
hist(w)


################################################################################
#                             Full data log-likelihood                         #
################################################################################

- full_data_log_likelihood(x = x, w = w, mu = oracle_param$mu + 1, sigma = oracle_param$sigma + 4, nu = oracle_param$nu)
- full_data_log_likelihood2(x = x, w = w, mu = oracle_param$mu, sigma = oracle_param$sigma, nu = oracle_param$nu)

mu_values <- seq(4.5, 6.5, length.out = 10)  # Adjust the range as needed
sigma_values <- seq(2, 4, length.out = 10)   # Avoid zero to prevent division by zero errors


FDLL_func_fac <- FDLL_function_factory(x = x, w = w, nu = oracle_param$nu)

FDLL_func_fac$H(c(oracle_param$mu + 1000, oracle_param$sigma + 1000))

# Create a grid for mu and sigma
grid <- expand.grid(mu = mu_values, sigma = sigma_values)

# Calculate marginal_log_likelihood for each combination of mu and sigma
grid$log_likelihood <- mapply(function(mu, sigma) {
  FDLL_func_fac$H(c(mu, sigma))
}, grid$mu, grid$sigma)


breaks_seq <- seq(min(grid$log_likelihood), max(grid$log_likelihood), length.out = 10)  # Adjust length.out for smaller intervals

# Alternatively, for a contour plot
ggplot(grid, aes(x = mu, y = sigma, z = log_likelihood)) +
  geom_contour_filled(breaks = breaks_seq) +
  geom_point(aes(x = oracle_param$mu, y = oracle_param$sigma), color = "red") +
  labs(title = "Contour Plot of Marginal Log-Likelihood", x = "Mu", y = "Sigma") +
  theme_bw()

GD_func_fac$H(par = c(5.6, 3.1))



full_data_mle(x = x, w = w, nu = oracle_param$nu)




################################################################################
#                             EM_algorithm.R tests                             #
################################################################################


# Test functions in EM_algorithm.R

alpha_prime <- alpha_prime_func(oracle_param$nu)
beta_prime <- beta_prime_func(x, oracle_param$nu, oracle_param$mu, oracle_param$sigma)

mu_hat_func(x = x, beta_prime = beta_prime)
sigma_hat_func(x = x, nu = oracle_param$nu, mu_hat = oracle_param$mu, alpha_prime = alpha_prime, beta_prime = beta_prime)










################################################################################
#                         Gradient_descent.R tests                             #
################################################################################


marginal_log_likelihood(x = x, mu = oracle_param$mu, 
                        sigma = oracle_param$sigma, 
                        nu = oracle_param$nu)



GD_func_fac <- GD_function_factory(x, nu)

mu_values <- seq(4.5, 6.5, length.out = 50)     # Adjust the range as needed
sigma_values <- seq(2, 4, length.out = 50)   # Avoid zero to prevent division by zero errors

# Create a grid for mu and sigma
grid <- expand.grid(mu = mu_values, sigma = sigma_values)

# Calculate marginal_log_likelihood for each combination of mu and sigma
grid$log_likelihood <- mapply(function(mu, sigma) {
  GD_func_fac$H(c(mu, sigma))
}, grid$mu, grid$sigma)


breaks_seq <- seq(min(grid$log_likelihood), max(grid$log_likelihood), length.out = 20)  # Adjust length.out for smaller intervals

# Alternatively, for a contour plot
ggplot(grid, aes(x = mu, y = sigma, z = log_likelihood)) +
  geom_contour_filled(breaks = breaks_seq) +
  geom_point(aes(x = oracle_param$mu, y = oracle_param$sigma), color = "red") +
  labs(title = "Contour Plot of Marginal Log-Likelihood", x = "Mu", y = "Sigma") +
  theme_bw()

GD_func_fac$H(par = c(5.6, 3.1))







oracle_param <- parameters(mu = mean_oracle, sigma = sigma_oracle, nu = nu)
oracle_par <- oracle_param$par
x <- simulate(oracle_param, N = N)$x
marginal_log_likelihood(x = x, mu = oracle_param$mu, 
                        sigma = oracle_param$sigma, 
                        nu = oracle_param$nu)
d_marginal_log_likelihood(x = x, mu = oracle_param$mu, 
                          sigma = oracle_param$sigma, 
                          nu = oracle_param$nu)


GD_func_fac <- GD_function_factory(x, nu)

GD_func_fac$grad(c(oracle_param$mu, 
                   sigma = oracle_param$sigma))
GD_func_fac$H(c(oracle_param$mu, 
                sigma = oracle_param$sigma))

# Test gradient descent
grad_desc(par = c(3,5),
          H = GD_func_fac$H, 
          grad = GD_func_fac$grad,
          clipping = TRUE)





################################################################################
#                         Fisher_information.R tests                           #
################################################################################
N <- 10000
nu <- 1
mean_oracle <- 2
sigma_oracle <- 4


oracle_param <- parameters(mu = mean_oracle, sigma = sigma_oracle, nu = nu)
oracle_par <- oracle_param$par
x <- simulate(oracle_param, N = N)$x
w <- simulate(oracle_param, N = N)$w


alpha_prime <- alpha_prime_func(oracle_param$nu)
mu_opt <- EM_algorithm(x, init_theta = c(2, 4), nu = nu, max_iter = 5000, tolerance = 1e-3)[1]
sigma_opt <- EM_algorithm(x, init_theta = c(2, 4), nu = nu, max_iter = 5000, tolerance = 1e-3)[2]


Qgrad(x = x, nu = nu, mu = mu_opt, sigma = sigma_opt, alpha_prime = alpha_prime, 
      beta_prime = beta_prime_func(x = x, nu = nu, mu = mu_opt, sigma = sigma_opt))
fisher_information_naive(mu = mu_opt, sigma = sigma_opt, nu = nu, alpha_prime = alpha_prime, x = x)
fisher_information_naive2(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu)
fisher_information_naive3(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu)
fisher_information(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu)

microbenchmark::microbenchmark(
  fisher_information_naive(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu, alpha_prime = alpha_prime),
  fisher_information_naive2(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu),
  fisher_information_naive3(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu),
  fisher_information(x = x, mu = mu_opt, sigma = sigma_opt, nu = nu)
)

mu_var <- (fisher_inf(mu = mu_opt, sigma = sigma_opt, nu = nu, alpha_prime = alpha_prime, x = x)/N)[1,1]
sigma_var <- (fisher_inf(mu = mu_opt, sigma = sigma_opt, nu = nu, alpha_prime = alpha_prime, x = x)/N)[2,2]

mu_opt + 1.96 * sqrt(mu_var) * c(-1, 1)
sigma_opt + 1.96 * sqrt(sigma_var) * c(-1, 1)








dQ_dmu <- function(mu, sigma, nu, alpha_prime, beta_prime, x){
  - alpha_prime / (nu * sigma^2) * sum((x - mu) / beta_prime)
}

dQ_dsigma <- function(mu, sigma, nu, alpha_prime, beta_prime, x){
  n <- length(x)
  - n/sigma + alpha_prime / (nu * sigma^3) * sum((x - mu)^2 / beta_prime)
}


-numDeriv::jacobian(d_marginal_log_likelihood2_par, x = c(mu_hat, sigma_hat), nu = nu, data = x)
fisher_information(x = x, mu = mu_hat, sigma = sigma_hat, nu = nu)
fisher_information_naive(x = x, par = c(mu_hat, sigma_hat), nu = nu, alpha_prime = alpha_prime)
fisher_information_naive2(x = x, mu = mu_hat, sigma = sigma_hat, nu = nu)
fisher_information_naive3(x = x, mu = mu_hat, sigma = sigma_hat, nu = nu)






