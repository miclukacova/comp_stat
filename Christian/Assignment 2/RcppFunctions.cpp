#include <Rcpp.h>
using namespace Rcpp;

// This is the Rcpp document for assignment 2 of the course "Computational Statistics" at the University of Copenhagen.
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector rcpp_target_distribution_pois(NumericVector y, NumericVector x, NumericVector z, NumericVector xz){

  const int n = y.size();
  const int m = x.size();
  
  NumericVector out(n);
  
  // y*xz - exp(y*x)
  
  double temp_sum = 0;
  
  for (int i = 0; i < n; ++i) {
    
    temp_sum = 0;
    
    for (int j = 0; j < m; ++j) {
      temp_sum += y[i] * xz[j] - exp(y[i] * x[j]);
    }
    out[i] = exp(temp_sum);
  }
  return out;
}




NumericVector rcpp_target_distribution_pois_double(double y, NumericVector x, NumericVector z, NumericVector xz){
  
  const int m = x.size();
  
  double temp_sum = 0;
  
  for (int j = 0; j < m; ++j) {
    temp_sum += y * xz[j] - exp(y * x[j]);
  }
  return temp_sum;
}



// Assuming the rcpp_target_distribution_pois function is already defined elsewhere
// NumericVector rcpp_target_distribution_pois(NumericVector y, NumericVector x, NumericVector z, NumericVector xz);

// [[Rcpp::export]]
List gaussianRcpp(int N, double mu, double sigma, double alpha_star, 
                        NumericVector x_poisson, NumericVector z_poisson, NumericVector xz_poisson) {
  
  // Initialize samples
  NumericVector samples(N);
  
  // Count number of rejections
  int rejection_count = 0;
  
  // Variables for the loop
  double x0, u;
  bool reject;
  NumericVector target_x0;
  double proposal_x0;
  
  // Loop over N samples
  for (int i = 0; i < N; ++i) {
    reject = true;
    
    while (reject) {
      // Sample from the proposal distribution (normal distribution)
      x0 = R::rnorm(mu, sigma);
      
      // Sample from uniform distribution
      u = R::runif(0.0, 1.0);
      
      // Evaluate target distribution at x0
      target_x0 = rcpp_target_distribution_pois(NumericVector::create(x0), x_poisson, z_poisson, xz_poisson);
      
      // Evaluate proposal distribution at x0 (normal PDF)
      proposal_x0 = R::dnorm(x0, mu, sigma, false);
      
      // Acceptance-rejection step
      reject = (u > alpha_star * target_x0[0] / proposal_x0);
      
      if (reject) {
        // Count rejections
        rejection_count += 1;
      }
    }
    
    // Store accepted sample
    samples[i] = x0;
  }
  
  // Estimate acceptance rate
  double alpha_hat = (double) N / (N + rejection_count);
  
  // Return the results as a list
  return List::create(
    Named("samples") = samples,
    Named("rejection_count") = rejection_count,
    Named("alpha_hat") = alpha_hat
  );
}










// NumericVector poisson_cpp(int N, double mu, double sigma, 
//                           NumericVector x, NumericVector z, NumericVector xz,
//                           double alpha_star) {
//   NumericVector samples(N);
//   
//   double x0;
//   bool reject;
//   
//   for(int i = 0; i < N; ++i) {
//     do {
//       x0 = R::rnorm(mu, sigma);
//       reject = R::runif(0, 1) > alpha_star * rcpp_target_distribution_pois_double(x0, x, z, xz)/
//         R::dnorm(x0, mu, sigma); // Acceptance-rejection step
//     } while(reject);
//     samples[i] = x0;
//   }
//   return samples;
// }


