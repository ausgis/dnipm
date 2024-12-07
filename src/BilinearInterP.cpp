#include <Rcpp.h>
#include <vector>
#include <algorithm>

// Function to compute the bilinear interpolation for a given x and y
double bilinearInterpolation(double x, double y,
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

  // Extract the four bounding points
  double x0 = xs(i - 1, 0);
  double x1 = xs(i, 0);
  double y0 = ys(0, j - 1);
  double y1 = ys(0, j);

  double z00 = zs(i - 1, j - 1);
  double z01 = zs(i - 1, j);
  double z10 = zs(i, j - 1);
  double z11 = zs(i, j);

  // Check for NA values
  bool has_na = false;
  if (Rcpp::NumericVector::is_na(z00) || Rcpp::NumericVector::is_na(z01) ||
      Rcpp::NumericVector::is_na(z10) || Rcpp::NumericVector::is_na(z11)) {
    has_na = true;
  }

  // If NA_rm is false and there is an NA value, return NA
  if (has_na && !NA_rm) {
    return NA_REAL;
  }

  // Compute the Lagrange basis polynomials
  double L0x = (x - x1) / (x0 - x1);
  double L1x = (x - x0) / (x1 - x0);
  double L0y = (y - y1) / (y0 - y1);
  double L1y = (y - y0) / (y1 - y0);

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
