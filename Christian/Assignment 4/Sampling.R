#Sampling

N <- 5000
omega <- 1 # Chaning this significantly changes the output. High values make the algorithm worse

true_par <- c(2, 5, 1, 2)
x_i <- exp(rnorm(N, mean = 0, sd = omega^2))
y_i <- log_logistic_dose_response_model(x_i, true_par) + rnorm(N, mean = 0, sd = 0.01)

rate <- decay_scheduler(gamma0 = 1, a = 2, gamma1 = 1e-2, n = 100)

SGD_tracer <- tracer(c("par", "k"), Delta = 0) 

init_par <- c(1,1,1,1)
sgd(par0 = init_par, 
    grad = gradient, 
    gamma = rate,
    N = N,
    epoch = batch,
    m = 500,
    maxiter = 100,
    sampler = sample,
    cb = SGD_tracer$tracer,
    x = x_i,
    y = y_i)



squared_error_single <- function(x, y, par){
  return(sum((y - f(x, par))^2))
}

squared_error <- function(x, y, alpha, beta, gamma, rho){
  param <- cbind(alpha, beta, gamma, rho)
  apply(param, 1, function(par) squared_error_single(x, y, par))
}

SGD_trace <- summary(SGD_tracer)
SGD_trace <- transform(
  SGD_trace,
  loss = squared_error(x_i, y_i, par.1, par.2, par.3, par.4),
  H_distance = abs(squared_error(x_i, y_i, true_par[1], true_par[2], true_par[3], true_par[4]) - squared_error(x_i, y_i, par.1, par.2, par.3, par.4))
)
tail(SGD_trace)


# Plot on log scale

ggplot(SGD_trace, aes(x = .time, y = H_distance)) +
  geom_line() +
  scale_y_log10() +
  labs(title = "Loss vs Time", x = "Time", y = "Loss")



#Test objects

SGD_tracer <- tracer(c("par", "n"), Delta = 0)
SGD_object <- SGD(par0 = init_par, 
                  grad = gradient, 
                  gamma = rate,
                  N = N,
                  epoch = adam(),
                  m = 500,
                  maxiter = 100,
                  sampler = sample,
                  cb = SGD_tracer$tracer,
                  x = x_i,
                  y = y_i)
SGD_object

summary(SGD_object)
print(SGD_object)




# Define an S3 class and its print method
my_class <- function(a) {
  structure(list(a = a), class = "my_class")
}

# Custom print method for my_class
print.my_class <- function(x) {
  cat("This is a custom object:\n")
  print(x$a)  # Print the 'a' element
}

# Create an object of class "my_class"
obj <- my_class(c(1, 2, 3))

# When evaluating obj, print.my_class() is called
obj



