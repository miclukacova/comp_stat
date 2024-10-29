#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector kernbinloopRcpp(int B, NumericVector x, double l, double delta) {
  
  int n = x.size();
  NumericVector w(B);
  
  for (int j = 0; j < n; ++j){
    int i = floor((x[j] - l) / delta + 0.5); //Removed + 1 as Cpp is 0 index based
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






// [[Rcpp::export]]
NumericVector cv_rcpp(NumericVector h_grid, NumericVector x, IntegerVector group, int folds) {
  int H = h_grid.size();
  int n = x.size();
  
  NumericVector ll_sum(H); // Vector to store the log-likelihood sums for each h
  
  for (int h_index = 0; h_index < H; h_index++) {
    double h = h_grid[h_index];
    NumericVector fold_l_cv(folds, 0.0); // Vector to store log-likelihood for each fold
    
    for (int fold = 1; fold <= folds; fold++) {
      std::vector<double> other_sets;
      std::vector<double> hold_out_set;
      
      // Define hold-out set and remaining sets for fold
      for (int i = 0; i < n; i++) {
        if (group[i] == fold) {
          hold_out_set.push_back(x[i]);
        } else {
          other_sets.push_back(x[i]);
        }
      }
      
      int N_i = other_sets.size();
      int M = hold_out_set.size();
      NumericMatrix kernel_input(M, N_i);
      
      // Precalculate kernel input matrix
      for (int j = 0; j < M; j++) {
        for (int k = 0; k < N_i; k++) {
          kernel_input(j, k) = (hold_out_set[j] - other_sets[k]) / h;
        }
      }
      
      // Apply the Epanechnikov kernel to the kernel input matrix
      NumericMatrix kernel_values = epanechnikovMatrixRcpp(kernel_input);
      
      // Calculate kernel density estimates for each hold-out point
      NumericVector f_i_hats(M, 0.0);
      for (int j = 0; j < M; j++) {
        for (int k = 0; k < N_i; k++) {
          f_i_hats[j] += kernel_values(j, k);
        }
        f_i_hats[j] /= (h * N_i);
      }
      
      // Calculate the log-likelihood for each fold
      for (int j = 0; j < M; j++) {
        fold_l_cv[fold - 1] += log(f_i_hats[j]);
      }
    }
    
    // Sum log-likelihoods to obtain log-likelihood for h
    ll_sum[h_index] = sum(fold_l_cv);
  }
  
  return ll_sum;
}





