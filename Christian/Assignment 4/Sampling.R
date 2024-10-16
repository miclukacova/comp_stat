SGD_tracer <- tracer(c("par", "k"), Delta = 0) 

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


#Sampling
set.seed(16102024)
N <- 10000
omega <- 0.5 # Chaning this significantly changes the output. High values make the algorithm worse

true_par <- c(2, 5, 1, 2)
init_par <- c(1,1,1,1)
x_i <- exp(rnorm(N, mean = 0, sd = omega^2))
y_i <- f(x_i, true_par) + rnorm(N, mean = 0, sd = 0.001)

iterations <- 1000
batch_size <- 20

# Decay schedule
rate_batch <- decay_scheduler(gamma0 = 1, a = 1, gamma1 = 1e-1, n = iterations)
rate_momentum <- decay_scheduler(gamma0 = 1, a = 1, gamma1 = 1e-1, n = iterations)
rate_adam <- decay_scheduler(gamma0 = 1e-1, a = 1, gamma1 = 1e-5, n = iterations)

#Test objects

# SGD_tracer_vanilla <- tracer(c("par", "n"), Delta = 0)
# SGD_object_vanilla <- SGD(par0 = init_par, grad = gradient, gamma = rate,
#                   N = N, epoch = NULL, maxiter = iterations, sampler = sample,
#                   cb = SGD_tracer_vanilla, x = x_i, y = y_i,
#                   true_par = true_par)

SGD_tracer_batch <- tracer(c("par", "n"), Delta = 0)
SGD_object_batch <- SGD(par0 = init_par, grad = gradient, gamma = rate_batch,
                       N = N, epoch = batch, m = batch_size, maxiter = iterations, sampler = sample,
                       cb = SGD_tracer_batch, x = x_i, y = y_i,
                       true_par = true_par)


SGD_tracer_momentum <- tracer(c("par", "n"), Delta = 0)
SGD_object_momentum <- SGD(par0 = init_par, grad = gradient, gamma = rate_momentum,
                        N = N, epoch = momentum(), m = batch_size, maxiter = iterations, sampler = sample,
                        cb = SGD_tracer_momentum, x = x_i, y = y_i,
                        true_par = true_par)


SGD_tracer_adam <- tracer(c("par", "n"), Delta = 0)
SGD_object_adam <- SGD(par0 = init_par, grad = gradient, gamma = rate_adam,
                       N = N, epoch = adam(), m = batch_size, maxiter = iterations, sampler = sample,
                       cb = SGD_tracer_adam, x = x_i, y = y_i,
                       true_par = true_par)

plot(SGD_object_adam, 3) + 
  geom_line(aes(x = plot_data(SGD_object_batch)$.time, 
                y = plot_data(SGD_object_batch)$abs_dist_from_par), col = "orange") + 
  # geom_line(aes(x = plot_data(SGD_object_vanilla)$.time, 
  #               y = plot_data(SGD_object_vanilla)$loss), col = "blue") + 
  geom_line(aes(x = plot_data(SGD_object_momentum)$.time, 
                y = plot_data(SGD_object_momentum)$abs_dist_from_par), col = "red") +
  xlim(0,5)







