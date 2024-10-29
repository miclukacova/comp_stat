########### Class ##############

my_density <- function(x) {
  structure(
    list(
      x = x,
      cv_bw = cv_bw_l_fast(x, do_parallel = TRUE, loop = inner_loop_rcpp, k = 20),
      amise_bw = AMISE_bw(x)
    ),
    class = "density_object"
  )
}

kern_dens <- function(x, h = "CV", binned = FALSE, B = 100) {
  UseMethod("kern_dens")
}

kern_dens.density_object <- function(object, h = "AMISE", binned = FALSE, B = 100) {
  if(h == "CV") h <- object$cv_bw
  if(h == "AMISE") h <- object$amise_bw
  if(binned) return(kern_dens_bin(x = object$x, h = h, B = B))
  kern_dens1(object$x, h)
}

plot.density_object <- function(obj, h, p = 1, binned = FALSE, B = 100) {
  dens_data <- as.data.frame(kern_dens(obj, h)[1:2])
  dens_r <- as.data.frame(density(obj$x, kernel = "epanechnikov", h)[1:2])
  
  if(binned == TRUE){
    dens_data_binned <- as.data.frame(kern_dens(obj, h, binned = TRUE, B = B)[1:2])
    if(p < 3) p <- 3
  }
  
  if(p == 1){
    pp <- ggplot(tibble(x = x), aes(x = obj$x)) + 
      geom_histogram(bins = 30, fill = "blue", alpha = 0.3, aes(x = obj$x, y = ..density..)) +
      geom_function(fun = x_dens, alpha = 0.5, aes(color = "True dens"), linewidth = 1.1)+
      geom_line(data = dens_data, aes(x = x, y = y, color = "Est dens"), linetype = "dashed", linewidth = 1.1)+
      scale_color_manual(values = c("True dens" = "blue", "Est dens" = "hotpink"))+
      labs(title = "Comparison of true and estimated density", x = "x", y = "Density", color = "Density types")
  }
  if(p == 2){
    pp <- ggplot(tibble(x = x), aes(x = obj$x)) + 
      geom_histogram(bins = 30, fill = "blue", alpha = 0.3, aes(x = obj$x, y = ..density..)) +
      geom_function(fun = x_dens, alpha = 0.5, aes(color = "True dens"), linewidth = 1.1)+
      geom_line(data = dens_data, aes(x = x, y = y, color = "Est dens"), linetype = "dashed", linewidth = 1.1)+
      geom_line(data = dens_r, aes(x = x, y = y, color = "density"), linetype = "dashed", linewidth = 1.1) +
      scale_color_manual(values = c("True dens" = "blue", "Est dens" = "hotpink", "density" = "green2"))+
      labs(title = "Comparison of true and estimated density", x = "x", y = "Density", color = "Density types")
  }
  
  if(p == 3){
    pp <- ggplot(tibble(x = x), aes(x = obj$x)) + 
      geom_histogram(bins = 30, fill = "blue", alpha = 0.3, aes(x = obj$x, y = ..density..)) +
      geom_function(fun = x_dens, alpha = 0.5, linewidth = 1.1, aes(color = "True dens"))+
      geom_line(data = dens_data, aes(x = x, y = y, color = "Est dens"), linetype = "dashed", linewidth = 1.1)+
      geom_line(data = dens_r, aes(x = x, y = y, color = "density"), linetype = "dashed", linewidth = 1.1) +
      geom_line(data = dens_data_binned, aes(x = x, y = y, color = "Binned"), linetype = "dashed", linewidth = 1.1) +
      scale_color_manual(values = c("True dens" = "blue", "Est dens" = "hotpink", 
                                    "density" = "green2", "Binned" = "purple"))+
      labs(title = "Comparison of true and estimated density", x = "x",
           y = "Density", color = "Density types")
  }
  
  if(p == 4){
    pp <- ggplot(tibble(x = x), aes(x = obj$x)) + 
      geom_histogram(bins = 30, fill = "blue", alpha = 0.3, aes(x = obj$x, y = ..density..)) +
      geom_function(fun = x_dens, alpha = 0.5, aes(color = "True dens"), linewidth = 1.1)+
      geom_line(data = dens_data_binned, aes(x = x, y = y, color = "Binned"), linetype = "dashed", linewidth = 1.1) +
      scale_color_manual(values = c("True dens" = "blue", "Binned" = "purple"))+
      labs(title = "Comparison of true and estimated density", x = "x",
           y = "Density", color = "Density types")
  }
  return(pp)
  
}

test_dens.density_object <- function(x) {
  UseMethod("test_dens")
}

test_dens <- function(obj, level = 1e-3, binned = FALSE){
  h <- seq(0.2,2, length.out = 4)
  
  for(i in seq_along(h)){
    
    if(binned == TRUE) dens_data <- as.data.frame(kern_dens(obj, h[i], binned = TRUE)[1:2])
    else dens_data <- as.data.frame(kern_dens(obj, h[i])[1:2])
    
    dens_r <- as.data.frame(density(obj$x, kernel = "epanechnikov", h[i])[1:2])
    
    test_that("Our kernel density estimate correspond to density", {
      expect_equal(
        dens_data$y,
        dens_r$y,
        tolerance = level)})
  }
}

plot_data <- function(obj, h, binned = FALSE, B = 100) {
  UseMethod("plot_data")
}

plot_data.density_object <- function(obj, h, binned = FALSE, B = 100){
  dens_data <- as.data.frame(kern_dens(obj, h, binned = binned, B)[1:2])
  dens_data
}

# Print method
print.density_object <- function(object){
  cat("The optimal CV band width:\n")
  print(object$cv_bw)
  cat("The optimal AMISE band width:\n")
  print(object$amise_bw)
  cat("Data length:\n")
  print(length(x))
}
