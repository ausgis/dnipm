// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppKNNIndice
Rcpp::IntegerVector RcppKNNIndice(double x, double y, Rcpp::NumericMatrix xys, int k);
RcppExport SEXP _dnipm_RcppKNNIndice(SEXP xSEXP, SEXP ySEXP, SEXP xysSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xys(xysSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppKNNIndice(x, y, xys, k));
    return rcpp_result_gen;
END_RCPP
}
// lagrangeInterp
Rcpp::NumericVector lagrangeInterp(Rcpp::NumericMatrix xy, Rcpp::NumericMatrix xys, Rcpp::NumericVector zs);
RcppExport SEXP _dnipm_lagrangeInterp(SEXP xySEXP, SEXP xysSEXP, SEXP zsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xys(xysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type zs(zsSEXP);
    rcpp_result_gen = Rcpp::wrap(lagrangeInterp(xy, xys, zs));
    return rcpp_result_gen;
END_RCPP
}
// bilinearInterp
Rcpp::NumericVector bilinearInterp(Rcpp::NumericMatrix xy, Rcpp::NumericMatrix xys, Rcpp::NumericVector zs);
RcppExport SEXP _dnipm_bilinearInterp(SEXP xySEXP, SEXP xysSEXP, SEXP zsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xys(xysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type zs(zsSEXP);
    rcpp_result_gen = Rcpp::wrap(bilinearInterp(xy, xys, zs));
    return rcpp_result_gen;
END_RCPP
}
// bicubicInterp
Rcpp::NumericVector bicubicInterp(Rcpp::NumericMatrix xy, Rcpp::NumericMatrix xys, Rcpp::NumericVector zs);
RcppExport SEXP _dnipm_bicubicInterp(SEXP xySEXP, SEXP xysSEXP, SEXP zsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xys(xysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type zs(zsSEXP);
    rcpp_result_gen = Rcpp::wrap(bicubicInterp(xy, xys, zs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dnipm_RcppKNNIndice", (DL_FUNC) &_dnipm_RcppKNNIndice, 4},
    {"_dnipm_lagrangeInterp", (DL_FUNC) &_dnipm_lagrangeInterp, 3},
    {"_dnipm_bilinearInterp", (DL_FUNC) &_dnipm_bilinearInterp, 3},
    {"_dnipm_bicubicInterp", (DL_FUNC) &_dnipm_bicubicInterp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dnipm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
