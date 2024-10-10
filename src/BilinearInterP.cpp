#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector BilinearInterP(Rcpp::NumericMatrix xy,
                                   Rcpp::NumericMatrix xys,
                                   Rcpp::NumericVector zs) {
  int n = xy.nrow(); // Number of interpolation points
  Rcpp::NumericVector result(n); // Vector to store results

  int nsize = xys.nrow();

  for (int k = 0; k < n; ++k) {
    double x_target = xy(k,1); // Target x value for interpolation
    double y_target = xy(k,2); // Target y value for interpolation

    // Find the grid cell containing
    int x1 = 0, x2 = 0, y1 = 0, y2 = 0;

    for (int i = 0; i < nsize - 1; ++i) {
      if (xys(i,1) <= x_target && xys(i+1,1) > x_target) {
        x1 = i; // Lower x index
        x2 = i + 1; // Upper x index
        break;
      }
    }

    for (int j = 0; j < nsize - 1; j++) {
      if (xys(j,1) <= y_target && xys(j,2) > y_target) {
        y1 = j; // Lower y index
        y2 = j + 1; // Upper y index
        break;
      }
    }

    // Perform bilinear interpolation
    double Q11 = z[y1 * ncol + x1]; // Value at (x1, y1)
    double Q12 = z[y2 * ncol + x1]; // Value at (x1, y2)
    double Q21 = z[y1 * ncol + x2]; // Value at (x2, y1)
    double Q22 = z[y2 * ncol + x2]; // Value at (x2, y2)

    double x1_x0 = x[x1]; // x coordinate of lower x index
    double x2_x0 = x[x2]; // x coordinate of upper x index
    double y1_y0 = y[y1]; // y coordinate of lower y index
    double y2_y0 = y[y2]; // y coordinate of upper y index

    // Calculate the interpolated value
    result[k] = (Q11 * (x2_x0 - x_target) * (y2_y0 - y_target) +
      Q21 * (x_target - x1_x0) * (y2_y0 - y_target) +
      Q12 * (x2_x0 - x_target) * (y_target - y1_y0) +
      Q22 * (x_target - x1_x0) * (y_target - y1_y0)) /
        ((x2_x0 - x1_x0) * (y2_y0 - y1_y0));
  }

  return result; // Return the interpolated values
}
