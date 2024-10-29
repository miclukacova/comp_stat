# The target densitity

poisson_data <- read_csv("poisson.csv")
x <- poisson_data$x
z <- poisson_data$z

# constant term
sum1 <- sum(x * z)

# f_star
f_star <- function(y) exp(y * sum1  - sum(exp(y * x)))

# f_star vectorized
f_star_vec <- Vectorize(f_star)

# f_star_prime
df1 <- function(y) exp(y * sum1  - sum(exp(y * x))) * (sum1  - sum(x * exp(y * x)))

# Ratio to minimize
ratio <- function(y, mu, sigma) 1 / (sqrt(2 * pi *sigma^2)) * 
  exp(- (y - mu)^2 / ( 2 * sigma^2) - y * sum1  + sum(exp(y * x)))

ratio_vec <- Vectorize(ratio)

# plot
plot_ratio <- function(y) ratio_vec(y, mu = mu_opt, sigma = 0.05)

