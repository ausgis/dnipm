#include <Rcpp.h>
#include <vector>
#include <cmath> // for std::pow

// Function to calculate the B-spline basis function B_{i,k}(x)
double bsplineBasis(double x, int i, int k,
                    const std::vector<double>& knots) {
  // Base case: degree 0 (piecewise constant)
  if (k == 0) {
    if (knots[i] <= x && x < knots[i + 1]) {
      return 1.0;
    } else {
      return 0.0;
    }
  }

  // Recursive case: degree k > 0
  double term1 = 0.0, term2 = 0.0;

  // Avoid division by zero
  if (knots[i + k] != knots[i]) {
    term1 = (x - knots[i]) / (knots[i + k] - knots[i]) * bsplineBasis(x, i, k - 1, knots);
  }

  if (knots[i + k + 1] != knots[i + 1]) {
    term2 = (knots[i + k + 1] - x) / (knots[i + k + 1] - knots[i + 1]) * bsplineBasis(x, i + 1, k - 1, knots);
  }

  return term1 + term2;
}

// [[Rcpp::export]]
double bsplineInterpolation(double x, double y,
                            const Rcpp::NumericMatrix& xs,
                            const Rcpp::NumericMatrix& ys,
                            const Rcpp::NumericMatrix& zs,
                            int polydegreeX, int polydegreeY,
                            bool NA_rm = true) {
  int n = xs.ncol(); // Number of columns in xs (i dimension)
  int m = ys.nrow(); // Number of rows in ys (j dimension)
  double result = 0.0;

  // Extract knots for x and y
  std::vector<double> knotsX(n + polydegreeX + 1);
  std::vector<double> knotsY(m + polydegreeY + 1);

  // Initialize knots for x (uniform spacing)
  for (size_t i = 0; i < knotsX.size(); ++i) { // Use size_t for vector size
    knotsX[i] = static_cast<double>(i) / (knotsX.size() - 1);
  }

  // Initialize knots for y (uniform spacing)
  for (size_t j = 0; j < knotsY.size(); ++j) { // Use size_t for vector size
    knotsY[j] = static_cast<double>(j) / (knotsY.size() - 1);
  }

  // Loop over i (columns of xs) and j (rows of ys)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      double zij = zs(j, i); // f(xi, yj) value at (j, i)

      // Skip if zij is NA and NA_rm is true
      if (NA_rm && Rcpp::NumericVector::is_na(zij)) {
        continue;
      }

      // Calculate B-spline basis functions for x and y
      double Bx = bsplineBasis(x, i, polydegreeX, knotsX);
      double By = bsplineBasis(y, j, polydegreeY, knotsY);

      // Accumulate the result
      result += Bx * By * zij;
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector bsplineInterp(Rcpp::NumericMatrix xy,
                                  Rcpp::NumericMatrix xs,
                                  Rcpp::NumericMatrix ys,
                                  Rcpp::NumericMatrix zs,
                                  int polydegreeX, int polydegreeY,
                                  bool NA_rm = true) {
  int n = xy.nrow();
  Rcpp::NumericVector z_pred(n);

  for (int p = 0; p < n; ++p) {
    double x = xy(p, 0);
    double y = xy(p, 1);
    z_pred[p] = bsplineInterpolation(x, y, xs, ys, zs,
                                     polydegreeX, polydegreeY, NA_rm);
  }

  return z_pred;
}
