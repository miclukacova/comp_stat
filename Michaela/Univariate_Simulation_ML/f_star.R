# The target densitity

poisson_data <- read.csv("~/R/comp_stat/Univariate_Simulation_ML/poisson.csv")
x <- poisson_data$x
z <- poisson_data$z

# constant term
sum1 <- sum(x * z)

# f_star
f_star <- function(y) exp(y * sum1  - sum(exp(y * x)))

# f_star_prime
df1 <- function(y) exp(y * sum1  - sum(exp(y * x))) * (sum1  - sum(x * exp(y * x)))

# f_star vectorized
f_star_vec <- Vectorize(f_star)
