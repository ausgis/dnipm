#include <Rcpp.h>

// Function to calculate binomial coefficient C(n, k)
double binomialCoefficient(int n, int k) {
  if (k > n - k) {
    k = n - k; // Take advantage of symmetry
  }
  double result = 1.0;
  for (int i = 0; i < k; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }
  return result;
}

// [[Rcpp::export]]
double bezierInterpolation(double x, double y,
                           Rcpp::NumericMatrix xys,
                           Rcpp::NumericVector zs,
                           bool NA_rm = true) {
  int n = xys.nrow(); // Number of rows in xys
  double result = 0.0;

  for (int i = 0; i < n; ++i) {
    // double xi = xys(i, 0); // x-coordinate of the i-th point
    // double yi = xys(i, 1); // y-coordinate of the i-th point
    double zi = zs[i];     // z-value at (xi, yi)

    // Skip if zi is NA and NA_rm is true
    if (NA_rm && Rcpp::NumericVector::is_na(zi)) {
      continue;
    }

    // Calculate the binomial coefficients
    double binom_x = binomialCoefficient(n - 1, i);
    double binom_y = binomialCoefficient(n - 1, i);

    // Calculate the polynomial terms
    double poly_x = binom_x * std::pow(1 - x, n - i - 1) * std::pow(x, i);
    double poly_y = binom_y * std::pow(1 - y, n - i - 1) * std::pow(y, i);

    // Accumulate the result
    result += poly_x * poly_y * zi;
  }

  return result;
}
