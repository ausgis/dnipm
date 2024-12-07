#include <Rcpp.h>
#include <vector>
#include <algorithm>

// Function to compute the Lagrange basis polynomial for a given x and y
double LagrangeBasisOneAxis(double x, const Rcpp::NumericVector& x_coords, int i) {
  int n = x_coords.size();
  double Li = 1.0;

  for (int k = 0; k < n; ++k) {
    if (k != i) {
      Li *= (x - x_coords[k]) / (x_coords[i] - x_coords[k]);
    }
  }

  return Li;
}

// Function to compute the bicubic interpolation for a given x and y
double bicubicInterpolation(double x, double y,
                            const Rcpp::NumericMatrix& xs,
                            const Rcpp::NumericMatrix& ys,
                            const Rcpp::NumericMatrix& zs,
                            bool NA_rm = true) {
  int n = xs.nrow();
  int m = xs.ncol();

  // Find the bounding indices for x and y
  int i = 0, j = 0;
  while (i < n && xs(i, 0) < x) ++i;
  while (j < m && ys(0, j) < y) ++j;

  // Adjust indices to ensure they are within bounds
  if (i == 0) i = 1;
  if (j == 0) j = 1;
  if (i == n) i = n - 1;
  if (j == m) j = m - 1;

  // Extract the 3x3 grid of points around (x, y)
  Rcpp::NumericVector x_coords = {xs(i - 1, 0), xs(i, 0), xs(i + 1, 0)};
  Rcpp::NumericVector y_coords = {ys(0, j - 1), ys(0, j), ys(0, j + 1)};

  Rcpp::NumericMatrix z_grid(3, 3);
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      z_grid(ii, jj) = zs(i - 1 + ii, j - 1 + jj);
    }
  }

  // Check for NA values
  bool has_na = false;
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      if (Rcpp::NumericVector::is_na(z_grid(ii, jj))) {
        has_na = true;
        break;
      }
    }
    if (has_na) break;
  }

  // If NA_rm is false and there is an NA value, return NA
  if (has_na && !NA_rm) {
    return NA_REAL;
  }

  // Compute the bicubic interpolation
  double z = 0.0;
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      if (!Rcpp::NumericVector::is_na(z_grid(ii, jj))) {
        double Li = LagrangeBasisOneAxis(x, x_coords, ii);
        double Lj = LagrangeBasisOneAxis(y, y_coords, jj);
        z += z_grid(ii, jj) * Li * Lj;
      }
    }
  }

  return z;
}

// [[Rcpp::export]]
Rcpp::NumericVector bicubicInterp(Rcpp::NumericMatrix xy,
                                  Rcpp::NumericMatrix xs,
                                  Rcpp::NumericMatrix ys,
                                  Rcpp::NumericMatrix zs,
                                  bool NA_rm = true) {
  int n = xy.nrow();
  Rcpp::NumericVector z_pred(n);

  for (int p = 0; p < n; ++p) {
    double x = xy(p, 0);
    double y = xy(p, 1);
    z_pred[p] = bicubicInterpolation(x, y, xs, ys, zs, NA_rm);
  }

  return z_pred;
}

// // Function to compute the bicubic interpolation for a given x and y
// double bicubicInterpolation(double x, double y,
//                             const Rcpp::NumericMatrix& xs,
//                             const Rcpp::NumericMatrix& ys,
//                             const Rcpp::NumericMatrix& zs) {
//   int n = xs.nrow();
//   int m = xs.ncol();
//
//   // Find the bounding indices for x and y
//   int i = 0, j = 0;
//   while (i < n && xs(i, 0) < x) ++i;
//   while (j < m && ys(0, j) < y) ++j;
//
//   // Adjust indices to ensure they are within bounds
//   if (i == 0) i = 1;
//   if (j == 0) j = 1;
//   if (i == n) i = n - 1;
//   if (j == m) j = m - 1;
//
//   // Extract the 3x3 grid of points around (x, y)
//   Rcpp::NumericVector x_coords = {xs(i - 1, 0), xs(i, 0), xs(i + 1, 0)};
//   Rcpp::NumericVector y_coords = {ys(0, j - 1), ys(0, j), ys(0, j + 1)};
//
//   Rcpp::NumericMatrix z_grid(3, 3);
//   for (int ii = 0; ii < 3; ++ii) {
//     for (int jj = 0; jj < 3; ++jj) {
//       z_grid(ii, jj) = zs(i - 1 + ii, j - 1 + jj);
//     }
//   }
//
//   // Compute the bicubic interpolation
//   double z = 0.0;
//   for (int ii = 0; ii < 3; ++ii) {
//     for (int jj = 0; jj < 3; ++jj) {
//       double Li = LagrangeBasisOneAxis(x, x_coords, ii);
//       double Lj = LagrangeBasisOneAxis(y, y_coords, jj);
//       z += z_grid(ii, jj) * Li * Lj;
//     }
//   }
//
//   return z;
// }
//
// // [[Rcpp::export]]
// Rcpp::NumericVector bicubicInterp(Rcpp::NumericMatrix xy,
//                                   Rcpp::NumericMatrix xs,
//                                   Rcpp::NumericMatrix ys,
//                                   Rcpp::NumericMatrix zs) {
//   int n = xy.nrow();
//   Rcpp::NumericVector z_pred(n);
//
//   for (int p = 0; p < n; ++p) {
//     double x = xy(p, 0);
//     double y = xy(p, 1);
//     z_pred[p] = bicubicInterpolation(x, y, xs, ys, zs);
//   }
//
//   return z_pred;
// }
