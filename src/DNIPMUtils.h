#ifndef DNIPMUtils_H
#define DNIPMUtils_H

#include <Rcpp.h>

Rcpp::IntegerVector RcppKNNIndice(double x, double y,
                                  Rcpp::NumericMatrix xys,
                                  int k);

bool RcppRowEqualMat(const Rcpp::NumericMatrix& mat1,
                     const Rcpp::NumericMatrix& mat2,
                     int row1, int row2);

#endif // DNIPMUtils_H
