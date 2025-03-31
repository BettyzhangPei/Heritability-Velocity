#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumS3(NumericVector par, NumericMatrix G, NumericVector n, NumericVector resid_y, NumericVector t) {
  // number of subjects
  int N = G.nrow();
  
  // initialize S3
  double S3 = 0;
  
  // Precompute cumulative sums
  IntegerVector cum_n(N + 1);
  cum_n[0] = 0;
  for (int i = 1; i <= N; ++i) {
    cum_n[i] = cum_n[i - 1] + n[i - 1];
  }
  
  // Loop through each pair of persons
  for (int i = 0; i < N - 1; ++i) {
    int index_i_start = cum_n[i];
    int index_i_end = cum_n[i + 1];
    
    for (int k = i + 1; k < N; ++k) {
      int index_k_start = cum_n[k];
      int index_k_end = cum_n[k + 1];
      
      double G_ik = G(i, k);
      double G_ik_par0 = G_ik * par[0];
      double G_ik_par1 = G_ik * par[1];
      
      // Loop through each pair of measurements
      for (int j = index_i_start; j < index_i_end; ++j) {
        for (int m = index_k_start; m < index_k_end; ++m) {
          
          double resid_product = resid_y[j] * resid_y[m];
          double t_product = t[j] * t[m];
          
          double diff = resid_product - (G_ik_par0 + G_ik_par1 * t_product);
          
          S3 += diff * diff;
        }
      }
    }
  }
  
  return 2 * S3;
}

