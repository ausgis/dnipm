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
                                   Rcpp::NumericVector zs) {
  int n = xy.nrow();
  int m = xys.nrow();
  Rcpp::NumericVector z_pred(n);

  for (int p = 0; p < n; ++p) {
    double x = xy(p, 0);
    double y = xy(p, 1);
    double z_sum = 0.0;

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        double Lij = lagrangeBasis(x, y, xys, i, j);
        z_sum += zs[i * m + j] * Lij;
      }
    }

    z_pred[p] = z_sum;
  }

  return z_pred;
}
