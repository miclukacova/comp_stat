#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tar_dens_cpp(NumericVector y, NumericVector x, double zx) {
  int n = y.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    double sum_exp_yx = 0;
    for (int j = 0; j < x.size(); j++) {
      sum_exp_yx += exp(y[i] * x[j]);  // Recalculate sum(exp(y[i] * x)) for each y[i]
    }
    result[i] = exp(y[i] * zx - sum_exp_yx);
  }
  
  return result;
}


// [[Rcpp::export]]
NumericVector tar_dens_log_difference_cpp(NumericVector y, NumericVector x, double zx) {
  int n = y.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    // Compute the sum of x * exp(y[i] * x)
    double sum_exp = 0;
    for (int j = 0; j < x.size(); j++) {
      sum_exp += x[j] * exp(y[i] * x[j]);
    }
    // Calculate the difference
    result[i] = zx - sum_exp;
  }
  
  return result;
}