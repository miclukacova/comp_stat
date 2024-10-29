#### S3 class

my_sample <- function(n, sampler) {
  sample <- sampler(n)
  structure(
    list(
      "Length" = n,
      "Sample" = sample$sample,
      "Accept" = sample$accept
    ),
    class = "sampling_object"
  )
}


plot.sampling_object <- function(obj){
  plot_data <- tibble(y = obj$Sample)

  ggplot(data = plot_data, ) + 
    geom_histogram(bins = 30, fill = "blue", alpha = 0.3, aes(x = y, y = ..density..)) +
    labs(title = "Histogram of sample", x = "y",
         y = "Density")
}

# Print method
print.sampling_object <- function(object){
  cat("Number of data points:\n")
  print(object$Length)
  cat("Acceptance rate:\n")
  print(object$Accept)
}
