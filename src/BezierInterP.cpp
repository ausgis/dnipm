#include <Rcpp.h>
#include <cmath> // for std::pow

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
                           const Rcpp::NumericMatrix& xs,
                           const Rcpp::NumericMatrix& ys,
                           const Rcpp::NumericMatrix& zs,
                           bool NA_rm = true) {
  int n = xs.ncol(); // Number of columns in xs (i dimension)
  int m = ys.nrow(); // Number of rows in ys (j dimension)
  double result = 0.0;

  // Loop over i (columns of xs) and j (rows of ys)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      double xi = xs(0, i); // x-coordinate at column i
      double yj = ys(j, 0); // y-coordinate at row j
      double zij = zs(j, i); // f(xi, yj) value at (j, i)

      // Skip if zij is NA and NA_rm is true
      if (NA_rm && Rcpp::NumericVector::is_na(zij)) {
        continue;
      }

      // Calculate the binomial coefficients
      double binom_i = binomialCoefficient(n - 1, i);
      double binom_j = binomialCoefficient(m - 1, j);

      // Calculate the polynomial terms
      double poly_x = binom_i * std::pow(1 - x, n - i - 1) * std::pow(x, i);
      double poly_y = binom_j * std::pow(1 - y, m - j - 1) * std::pow(y, j);

      // Accumulate the result
      result += poly_x * poly_y * zij;
    }
  }

  return result;
}
