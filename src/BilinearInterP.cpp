#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <limits>

// Function to compute the bilinear interpolation for a given x and y
double bilinearInterpolation(double x, double y,
                             const Rcpp::NumericMatrix& xs,
                             const Rcpp::NumericMatrix& ys,
                             const Rcpp::NumericMatrix& zs,
                             bool NA_rm = true) {
  // Get dimensions of the input matrices
  int n = xs.nrow(); // Number of rows in xs
  int m = xs.ncol(); // Number of columns in xs

  // Check if x or y is out of bounds
  if (x < xs(0, 0) || x > xs(0, m - 1) || y < ys(0, 0) || y > ys(n - 1, 0)) {
    return std::numeric_limits<double>::quiet_NaN(); // Return NaN if out of bounds
  }

  // Find the bounding indices for x and y
  int i = 0, j = 0;
  while (i < n && ys(i, 0) < y) ++i; // Find the row index for y
  while (j < m && xs(0, j) < x) ++j; // Find the column index for x

  // Adjust indices to ensure they are within bounds
  if (i == 0) i = 1; // Ensure i is at least 1
  if (j == 0) j = 1; // Ensure j is at least 1
  if (i == n) i = n - 1; // Ensure i is at most n-1
  if (j == m) j = m - 1; // Ensure j is at most m-1

  // Extract the four bounding points
  double x0 = xs(0, j - 1); // Lower x bound
  double x1 = xs(0, j);     // Upper x bound
  double y0 = ys(i - 1, 0); // Lower y bound
  double y1 = ys(i, 0);     // Upper y bound

  // Extract the four z values at the bounding points
  double z00 = zs(i - 1, j - 1);
  double z01 = zs(i - 1, j);
  double z10 = zs(i, j - 1);
  double z11 = zs(i, j);

  // Check for NA values in the z values
  bool has_na = false;
  if (Rcpp::NumericVector::is_na(z00) || Rcpp::NumericVector::is_na(z01) ||
      Rcpp::NumericVector::is_na(z10) || Rcpp::NumericVector::is_na(z11)) {
    has_na = true;
  }

  // If NA_rm is false and there is an NA value, return NaN
  if (has_na && !NA_rm) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Compute the Lagrange basis polynomials
  double L0x = (x - x1) / (x0 - x1); // Lagrange basis for x0
  double L1x = (x - x0) / (x1 - x0); // Lagrange basis for x1
  double L0y = (y - y1) / (y0 - y1); // Lagrange basis for y0
  double L1y = (y - y0) / (y1 - y0); // Lagrange basis for y1

  // Compute the interpolated value
  double z = 0.0;
  if (!Rcpp::NumericVector::is_na(z00)) z += z00 * L0x * L0y;
  if (!Rcpp::NumericVector::is_na(z01)) z += z01 * L0x * L1y;
  if (!Rcpp::NumericVector::is_na(z10)) z += z10 * L1x * L0y;
  if (!Rcpp::NumericVector::is_na(z11)) z += z11 * L1x * L1y;

  return z;
}

// [[Rcpp::export]]
Rcpp::NumericVector bilinearInterp(Rcpp::NumericMatrix xy,
                                   Rcpp::NumericMatrix xs,
                                   Rcpp::NumericMatrix ys,
                                   Rcpp::NumericMatrix zs,
                                   bool NA_rm = true) {
  int n = xy.nrow();
  Rcpp::NumericVector z_pred(n);

  for (int p = 0; p < n; ++p) {
    double x = xy(p, 0);
    double y = xy(p, 1);
    z_pred[p] = bilinearInterpolation(x, y, xs, ys, zs, NA_rm);
  }

  return z_pred;
}

// // Function to compute the bilinear interpolation for a given x and y
// double bilinearInterpolation(double x, double y,
//                              const Rcpp::NumericMatrix& xs,
//                              const Rcpp::NumericMatrix& ys,
//                              const Rcpp::NumericMatrix& zs) {
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
//   // Extract the four bounding points
//   double x0 = xs(i - 1, 0);
//   double x1 = xs(i, 0);
//   double y0 = ys(0, j - 1);
//   double y1 = ys(0, j);
//
//   double z00 = zs(i - 1, j - 1);
//   double z01 = zs(i - 1, j);
//   double z10 = zs(i, j - 1);
//   double z11 = zs(i, j);
//
//   // Compute the Lagrange basis polynomials
//   double L0x = (x - x1) / (x0 - x1);
//   double L1x = (x - x0) / (x1 - x0);
//   double L0y = (y - y1) / (y0 - y1);
//   double L1y = (y - y0) / (y1 - y0);
//
//   // Compute the interpolated value
//   double z = z00 * L0x * L0y +
//     z01 * L0x * L1y +
//     z10 * L1x * L0y +
//     z11 * L1x * L1y;
//
//   return z;
// }
//
// // [[Rcpp::export]]
// Rcpp::NumericVector bilinearInterp(Rcpp::NumericMatrix xy,
//                                    Rcpp::NumericMatrix xs,
//                                    Rcpp::NumericMatrix ys,
//                                    Rcpp::NumericMatrix zs) {
//   int n = xy.nrow();
//   Rcpp::NumericVector z_pred(n);
//
//   for (int p = 0; p < n; ++p) {
//     double x = xy(p, 0);
//     double y = xy(p, 1);
//     z_pred[p] = bilinearInterpolation(x, y, xs, ys, zs);
//   }
//
//   return z_pred;
// }
