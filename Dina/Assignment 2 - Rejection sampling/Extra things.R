This also becomes clear if you plot it.
```{r, echo=FALSE, fig.align='center', fig.width= 7.5, fig.height=4}
# Assuming poisson_data is already defined and contains columns 'x' and 'z'
zx <- sum(poisson_data$x * poisson_data$z)

# Define the log of the target density function
tar_dens_log <- function(y) {
  sapply(y, function(yi) log(tar_dens(yi)))
}

# Generate a sequence of y values
y_vals <- seq(0, 1, length.out = 100)
log_dens_vals <- tar_dens_log(y_vals)

# Create a data frame for ggplot
plot_data <- data.frame(y = y_vals, log_density = log_dens_vals)

# Create the ggplot
ggplot(plot_data, aes(x = y, y = log_density)) +
  geom_line(color = base, size = 1) +  # Line for the log target density
  labs(
    x = "y",
    y = "log(f(y))"
  ) +
  my_theme

```