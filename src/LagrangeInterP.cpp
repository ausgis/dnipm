#include <Rcpp.h>
#include <vector>
#include <algorithm>

// Function to compute the Lagrange basis polynomial for a given x and y
// [[Rcpp::export]]
double lagrangeBasis(double x, double y,
                     const Rcpp::NumericMatrix& xys,
                     int i) {
  double Lij = 1.0; // Initialize the Lagrange basis value
  int m = xys.nrow(); // Number of known points

  for (int j = 0; j < m; ++j) {
    if (j != i) {
      double xi = xys(i, 0); // x-coordinate of the i-th known point
      double yi = xys(i, 1); // y-coordinate of the i-th known point
      double xj = xys(j, 0); // x-coordinate of the j-th known point
      double yj = xys(j, 1); // y-coordinate of the j-th known point

      // Check if denominator is zero
      double denom_x = xi - xj;
      double denom_y = yi - yj;

      if (denom_x == 0 || denom_y == 0) {
        // If denominator is zero, skip this iteration
        continue;
      }

      // Compute the Lagrange basis
      Lij *= (x - xj) / denom_x; // * (y - yj) / denom_y;
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
    for (int j = 0; j < m; ++j) {
      double Lij = lagrangeBasis(x, y, xys, j); // Compute Lagrange basis
      Rcpp::Rcout << "The value of Lij : " << Lij << "\n";
      double z_val = zs[j]; // Known z-value at (j, j)

      // If z_val is NA and NA_rm is true, skip this point
      if (Rcpp::NumericVector::is_na(z_val) && NA_rm) {
        continue;
      } else {
        z_sum += z_val * Lij; // Add to the interpolation sum
        valid_interpolation = true; // Mark interpolation as valid
      }
    }

    // If interpolation is valid, assign z_sum to z_pred[p]; otherwise, keep it as NA
    if (valid_interpolation) {
      z_pred[p] = z_sum;
    }
  }

  return z_pred;
}

// // Function to compute the Lagrange basis polynomial for a given x and y
// // [[Rcpp::export]]
// double lagrangeBasis(double x, double y,
//                      const Rcpp::NumericMatrix& xys,
//                      int j) {
//   double Lij = 1.0; // Initialize the Lagrange basis value
//   int m = xys.nrow(); // Number of known points
//   bool valid_multiplication = false; // Flag to track if any valid multiplication has occurred
//
//   for (int k = 0; k < m; ++k) {
//     if (k != j) {
//       double xk = xys(k, 0); // x-coordinate of the k-th known point
//       double yk = xys(k, 1); // y-coordinate of the k-th known point
//
//       // Check if denominator is zero
//       double denom_x = xys(j, 0) - xk; // Denominator for x
//       double denom_y = xys(j, 1) - yk; // Denominator for y
//
//       if (denom_x == 0 || denom_y == 0) {
//         // If denominator is zero, skip this iteration
//         continue;
//       }
//
//       // Perform the multiplication if denominators are valid
//       Lij *= (x - xk) / denom_x * (y - yk) / denom_y;
//       valid_multiplication = true; // Mark that a valid multiplication has occurred
//     }
//   }
//
//   // If no valid multiplication occurred (e.g., all denominators were zero), set Lij to 0.0
//   if (!valid_multiplication) {
//     Lij = 0.0;
//   }
//
//   return Lij;
// }
