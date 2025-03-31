#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SumD_c(NumericMatrix G, NumericVector resid_y, NumericVector n, NumericVector t) {
  // number of subjects
  int N = G.nrow();
  
  // initialize S3
  NumericVector S3(5);
  
  // Precompute the cumulative sums
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

      double G_ik_sq = G(i, k) * G(i, k);
      double G_ik = G(i, k);

      // Loop through each pair of measurements
      for (int j = index_i_start; j < index_i_end; ++j) {
        for (int m = index_k_start; m < index_k_end; ++m) {
          
          double t_product = t[j] * t[m];
          double t_product_squared = t_product * t_product;
          double resid_product = resid_y[j] * resid_y[m];
          double resid_t_product = resid_product * t_product;

          S3[0] += G_ik_sq;
          S3[1] += G_ik_sq * t_product;
          S3[2] += G_ik_sq * t_product_squared;
          S3[3] += G_ik * resid_product;
          S3[4] += G_ik * resid_t_product;
        }
      }
    }
  }

  return 2 * S3;
}

