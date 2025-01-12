#include <Rcpp.h>
#include <vector>
#include <algorithm>

// Function to compute the Lagrange basis polynomial for a given x and y
double lagrangeBasis(double x, double y,
                     const Rcpp::NumericMatrix& xys,
                     int i, int j) {
  int n = xys.nrow();
  double Lij = 1.0;

  for (int k = 0; k < n; ++k) {
    if (k != i) {
      Lij *= (x - xys(k, 0)) / (xys(i, 0) - xys(k, 0));
    }
  }

  for (int k = 0; k < n; ++k) {
    if (k != j) {
      Lij *= (y - xys(k, 1)) / (xys(j, 1) - xys(k, 1));
    }
  }

  return Lij;
}

// [[Rcpp::export]]
Rcpp::NumericVector lagrangeInterp(Rcpp::NumericMatrix xy,
                                   Rcpp::NumericMatrix xys,
                                   Rcpp::NumericVector zs,
                                   bool NA_rm = true) {
  int n = xy.nrow(); // Number of points to interpolate
  int m = xys.nrow(); // Number of known points
  Rcpp::NumericVector z_pred(n, Rcpp::NumericVector::get_na()); // Initialize z_pred with NA values

  // Loop through each point to interpolate
  for (int p = 0; p < n; ++p) {
    double x = xy(p, 0); // x-coordinate of the point to interpolate
    double y = xy(p, 1); // y-coordinate of the point to interpolate
    double z_sum = 0.0; // Sum of interpolated values
    bool valid_interpolation = false; // Flag to track if interpolation is valid

    // Loop through known points to compute Lagrange basis and interpolate
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        double Lij = lagrangeBasis(x, y, xys, i, j); // Compute Lagrange basis
        double z_val = zs[i * m + j]; // Known z-value at (i, j)

        // If z_val is NA and NA_rm is true, skip this point
        if (Rcpp::NumericVector::is_na(z_val) && NA_rm) {
          continue;
        } else {
          z_sum += z_val * Lij; // Add to the interpolation sum
          valid_interpolation = true; // Mark interpolation as valid
        }
      }
    }

    // If interpolation is valid, assign z_sum to z_pred[p]; otherwise, keep it as NA
    if (valid_interpolation) {
      z_pred[p] = z_sum;
    }
  }

  return z_pred;
}

// // [[Rcpp::export]]
// Rcpp::NumericVector lagrangeInterp(Rcpp::NumericMatrix xy,
//                                    Rcpp::NumericMatrix xys,
//                                    Rcpp::NumericVector zs,
//                                    bool NA_rm = true) {
//   int n = xy.nrow();
//   int m = xys.nrow();
//   Rcpp::NumericVector z_pred(n);
//
//   for (int p = 0; p < n; ++p) {
//     double x = xy(p, 0);
//     double y = xy(p, 1);
//     double z_sum = 0.0;
//
//     for (int i = 0; i < m; ++i) {
//       for (int j = 0; j < m; ++j) {
//         double Lij = lagrangeBasis(x, y, xys, i, j);
//         double z_val = zs[i * m + j];
//
//         if (Rcpp::NumericVector::is_na(z_val) && NA_rm) {
//           break;
//         } else {
//           z_sum += z_val * Lij;
//         }
//       }
//     }
//
//     z_pred[p] = z_sum;
//   }
//
//   return z_pred;
// }
