# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

bicubicInterp <- function(xy, xs, ys, zs, NA_rm = TRUE) {
    .Call('_dnipm_bicubicInterp', PACKAGE = 'dnipm', xy, xs, ys, zs, NA_rm)
}

bilinearInterp <- function(xy, xs, ys, zs, NA_rm = TRUE) {
    .Call('_dnipm_bilinearInterp', PACKAGE = 'dnipm', xy, xs, ys, zs, NA_rm)
}

lagrangeInterp <- function(xy, xys, zs, NA_rm = TRUE) {
    .Call('_dnipm_lagrangeInterp', PACKAGE = 'dnipm', xy, xys, zs, NA_rm)
}

