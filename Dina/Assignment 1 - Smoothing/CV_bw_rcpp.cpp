#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Epanechnikov kernel function
arma::mat kernel(const arma::mat& z) {
  arma::mat abs_z = arma::abs(z);
  arma::mat result = arma::zeros(z.n_rows, z.n_cols);
  // Apply the Epanechnikov kernel formula where abs(z) <= 1
  result.elem(arma::find(abs_z <= 1)) = 0.75 * (1 - arma::square(z.elem(arma::find(abs_z <= 1))));
  return result;
}

// [[Rcpp::export]]
arma::vec calculate_log_likelihood(const arma::vec& x, const arma::vec& group, 
                                   const arma::vec& h, int k) {
  int n_h = h.n_elem;
  arma::vec ll(n_h, arma::fill::zeros);  // log-likelihood values for each h
  
  // Precompute inverses of h values
  arma::vec inv_h = 1 / h;
  
  // Loop over each bandwidth value in h
  for (int j = 0; j < n_h; ++j) {
    arma::vec f_hat;  // Initialize empty f_hat vector
    
    // Loop over each fold
    for (int i = 1; i <= k; ++i) {
      arma::uvec idx_in_group = arma::find(group != i);    // Indices not in group i
      arma::uvec idx_out_group = arma::find(group == i);   // Indices in group i
      
      int N_i = idx_in_group.n_elem;  // Number of elements not in group i
      arma::vec x_in_group = x.elem(idx_in_group);
      arma::vec x_out_group = x.elem(idx_out_group);
      
      // Optimize x_matrix calculation using precomputed inv_h and broadcasting
      arma::mat x_matrix = inv_h[j] * (x_in_group * arma::ones<arma::rowvec>(x_out_group.n_elem) - 
        arma::ones<arma::vec>(N_i) * x_out_group.t());
      
      // Evaluate the Epanechnikov kernel at each point
      arma::mat kerns = kernel(x_matrix);
      
      // Sum kernel values over columns
      arma::vec kern_sums = arma::sum(kerns, 0).t(); // Sum over columns
      
      // Calculate f_hat_i for this group and store in f_hat
      arma::vec f_hat_i = (1.0 / (h[j] * N_i)) * kern_sums;
      f_hat = arma::join_cols(f_hat, f_hat_i);
    }
    
    // Calculate log-likelihood for current h
    ll[j] = arma::sum(arma::log(f_hat));
  }
  
  return ll;  // Return log-likelihood vector
}
