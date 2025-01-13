#ifndef LAGRANGEBASIS_H
#define LAGRANGEBASIS_H

#include <Rcpp.h>

double lagrangeBasis(double x, const Rcpp::NumericVector& x_coords, int i);

#endif // LAGRANGEBASIS_H
