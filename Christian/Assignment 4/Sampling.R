#Sampling

N <- 20000
omega <- 1

x_i <- exp(rnorm(N, mean = 0, sd = omega^2))
y_i <- log_logistic_dose_response_model(x, c(5, 3, 1, 2)) + rnorm(N, mean = 0, sd = 1)

init_par <- c(1,1,1,1)
sgd(par = init_par, 
    grad = gradient, 
    gamma = 0.1,
    N = N,
    epoch = adam(),
    maxiter = 1000,
    sampler = sample,
    cb = NULL,
    x = x_i,
    y = y_i)
