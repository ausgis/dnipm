#include <Rcpp.h>
#include "DNIPMUtils.h"
using namespace Rcpp;

Rcpp::NumericVector LagrangeBasis(double x,
                                  Rcpp::NumericVector xs) {
  Rcpp::NumericVector res (xs.size(),1);
  for(int i = 0; i < xs.size(); ++i){
    double basis = 1.0;
    for (int j = 0; j < xs.size(); ++j) {
      if (j != i) {
        basis *= (x - xs[j]) / (xs[i] - xs[j]);
      }
    }
    res[i] = basis;
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector LagrangeInterp(Rcpp::NumericMatrix xy,
                                   Rcpp::NumericMatrix xys,
                                   Rcpp::NumericVector zs) {
  Rcpp::NumericVector res (xy.nrow());

  for (int n = 0; n < xy.nrow(); ++n){
    double z_interp = 0.0;
    Rcpp::NumericVector Lx = LagrangeBasis(xy(n,1),xys(_,1));
    Rcpp::NumericVector Ly = LagrangeBasis(xy(n,2),xys(_,2));
    for(int i = 0; i < xys.nrow(); ++i){
      z_interp += zs[i] * Lx[i] * Ly[i];
    }
    res[n] = z_interp;
  }
  return res;
}
