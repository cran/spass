// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// minFunc
double minFunc(NumericVector x, NumericVector daten, int dataNA);
RcppExport SEXP spass_minFunc(SEXP xSEXP, SEXP datenSEXP, SEXP dataNASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type daten(datenSEXP);
    Rcpp::traits::input_parameter< int >::type dataNA(dataNASEXP);
    __result = Rcpp::wrap(minFunc(x, daten, dataNA));
    return __result;
END_RCPP
}
// minFuncMult
double minFuncMult(NumericVector x, NumericMatrix daten, NumericVector dataNA, int n);
RcppExport SEXP spass_minFuncMult(SEXP xSEXP, SEXP datenSEXP, SEXP dataNASEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type daten(datenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dataNA(dataNASEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(minFuncMult(x, daten, dataNA, n));
    return __result;
END_RCPP
}
// minFuncBlinded
double minFuncBlinded(NumericVector x, NumericMatrix daten, NumericVector dataNA, NumericVector n, double delta);
RcppExport SEXP spass_minFuncBlinded(SEXP xSEXP, SEXP datenSEXP, SEXP dataNASEXP, SEXP nSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type daten(datenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dataNA(dataNASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    __result = Rcpp::wrap(minFuncBlinded(x, daten, dataNA, n, delta));
    return __result;
END_RCPP
}
