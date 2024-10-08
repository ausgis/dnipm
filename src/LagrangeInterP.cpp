#include <Rcpp.h>
using namespace Rcpp;

double LagrangeBasis(double x, int i, NumericVector xs) {
  double basis = 1.0;
  int pn = xs.size();

  // Loop through all points to calculate the basis polynomial
  for (int j = 0; j < pn; ++j) {
    if (j != i) {
      basis *= (x - xs[j]) / (xs[i] - xs[j]);
    }
  }
  return basis;
}

// [[Rcpp::export]]
double LagrangeInterp(double x, double y,
                      Rcpp::NumericMatrix xy,
                      Rcpp::NumericVector zs) {
  double res = 0.0;
  int n = xy.ncol();

  for(int i = 0; i < n; ++i){
    double Lx = LagrangeBasis(x, i ,xy(_,1));
    double Ly = LagrangeBasis(y, i ,xy(_,2));
    res += Lx * Ly * zs[i];
  }
  return res;
}
