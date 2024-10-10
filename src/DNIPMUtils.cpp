#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector RcppKNNIndice(double x, double y, Rcpp::NumericMatrix xys, int k) {
  int n = xys.nrow();
  int num_neighbors = (k >= 3) ? (k * k - 1) : 4;
  std::vector<std::pair<double, int>> distances(n);

  // Calculate the distance from each point in xys to the point (x, y)
  for (int i = 0; i < n; ++i) {
    double dx = xys(i, 0) - x;
    double dy = xys(i, 1) - y;
    double dist = std::sqrt(dx * dx + dy * dy);
    distances[i] = std::make_pair(dist, i);
  }

  // Sort the distances
  std::sort(distances.begin(), distances.end());

  // Select the indices of the closest num_neighbors points
  Rcpp::IntegerVector indices(num_neighbors);
  for (int i = 0; i < num_neighbors; ++i) {
    indices[i] = distances[i].second;
  }

  return indices;
}

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
