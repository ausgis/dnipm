#' two-dimensional numerical interpolation
#'
#' @param z The attribute vector
#' @param coords The coordinate matrix
#' @param predicts Matrix of locations to predict
#' @param method (optional) Interpolation method
#' @param polyorder (optional) Degree of the polynomial
#' @param na.rm (optional) Whether `NA` values should be stripped
#'
#' @return A numeric vector
#' @export
#'
#' @examples
#' data("ndvi")
#' interp(ndvi$NDVI,ndvi[,c("X","Y")],
#'        matrix(c(2745,2455),ncol = 2))
#'
interp = \(z, coords, predicts, method = "bilinear",
           polyorder = 3, na.rm = TRUE){
  predicts = as.matrix(predicts)
  coords = as.matrix(coords)
  xyzs = .xyz_vec2mat(z,coords)
  if (method == "lagrange") {
    res = lagrangeInterp(predicts,coords,z,na.rm)
  } else if (method == "bilinear") {
    res = bilinearInterp(predicts,xyzs[[2]],xyzs[[3]],xyzs[[1]],na.rm)
  } else if (method == "cubic") {
    res = bicubicInterp(predicts,xyzs[[2]],xyzs[[3]],xyzs[[1]],na.rm)
  } else if (method == "bezier") {
    res = bezierInterp(predicts,xyzs[[2]],xyzs[[3]],xyzs[[1]],na.rm)
  }else if (method == "bspline") {
    res = bsplineInterp(predicts,xyzs[[2]],xyzs[[3]],xyzs[[1]],
                        polyorder, polyorder, na.rm)
  }
  return(res)
}
