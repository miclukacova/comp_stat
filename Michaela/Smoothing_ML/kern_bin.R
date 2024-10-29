#-------------------------------------------------------------------------------
# Function computing the percentage of observations in each bin
#-------------------------------------------------------------------------------

kern_bin <- function(x, l, u, B) {
  w <- numeric(B)
  delta <- (u - l) / (B - 1)
  for (j in seq_along(x)) {
    i <- floor((x[j] - l) / delta + 0.5) + 1
    w[i] <- w[i] + 1
  }
  w / sum(w)
}


kern_bin_fast <- function(x, l, u, B) {
  delta <- (u - l) / (B - 1)
  is <- floor((x - l) / delta + 0.5) + 1
  w <- tabulate(is)
  w/sum(w)
}