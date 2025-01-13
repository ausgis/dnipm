#include <Rcpp.h>

// Function to compute the Lagrange basis polynomial for a given x
// [[Rcpp::export]]
double lagrangeBasis(double x, const Rcpp::NumericVector& x_coords, int i) {
  int n = x_coords.size();
  double Li = 1.0;

  for (int k = 0; k < n; ++k) {
    if (k != i) {
      // Check if denominator is zero
      double denom_x = x_coords[i] - x_coords[k];

      if (denom_x == 0) {
        // If denominator is zero, skip this iteration
        continue;
      }

      Li *= (x - x_coords[k]) / (denom_x);
    }
  }

  return Li;
}
