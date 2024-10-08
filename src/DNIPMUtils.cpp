#include <Rcpp.h>
using namespace Rcpp;

bool RcppRowEqualMat(const Rcpp::NumericMatrix& mat1,
                     const Rcpp::NumericMatrix& mat2,
                     int row1, int row2) {
  // Get the number of columns for both matrices
  int n_cols1 = mat1.ncol();
  int n_cols2 = mat2.ncol();

  // Check if the number of columns is the same in both matrices
  if (n_cols1 != n_cols2) {
    Rcpp::stop("Number of columns in the matrices do not match");
  }

  // Compare the elements of row1 from mat1 and row2 from mat2
  for (int j = 0; j < n_cols1; ++j) {
    if (mat1(row1, j) != mat2(row2, j)) {
      return false; // Return false if any element is different
    }
  }

  // Return true if all elements in the rows are the same
  return true;
}
