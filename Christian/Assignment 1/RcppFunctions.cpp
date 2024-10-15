#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector kernbinloopRcpp(int B, NumericVector x, double l, double delta) {
  
  int n = x.size();
  NumericVector w(B);
  
  for (int j = 0; j < n; ++j){
    int i = floor((x[j] - l) / delta + 0.5); //He removed + 1 as Cpp is 0 index based
    w[i] = w[i] + 1;
  }
  return w;
}



// [[Rcpp::export]]
double matrix_sum_pmaxRcpp(double r_hat, NumericMatrix abs_val_matrix) {
  double sum_val = 0.0;
  
  // Loop over the matrix and apply the formula
  for (int i = 0; i < abs_val_matrix.nrow(); i++) {
    for (int j = 0; j < abs_val_matrix.ncol(); j++) {
      // Compute the pmax(0, 2 * r_hat + abs_val_matrix[i,j])
      sum_val += std::max(0.0, 2 * r_hat + abs_val_matrix(i, j));
    }
  }
  
  return sum_val;
}

// [[Rcpp::export]]
NumericMatrix epanechnikovMatrixRcpp(NumericMatrix x){
  int n = x.nrow();
  int m = x.ncol();
  NumericMatrix out(n, m);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      double val = x(i, j);
      
      // Check the condition abs(x) <= 1
      if (std::abs(val) <= 1) {
        // Apply the formula (3/4) * (1 - x^2)
        out(i, j) = (3.0 / 4.0) * (1 - val * val);
      } else {
        // Set to 0 if the condition abs(x) > 1
        out(i, j) = 0.0;
      }
    }
  }
    return out;
}



