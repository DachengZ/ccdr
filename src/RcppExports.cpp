// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gridCCDr
List gridCCDr(NumericVector cors, List init_betas, IntegerVector nj, NumericVector aj, NumericVector lambdas, NumericVector params, int verbose);
RcppExport SEXP ccdr_gridCCDr(SEXP corsSEXP, SEXP init_betasSEXP, SEXP njSEXP, SEXP ajSEXP, SEXP lambdasSEXP, SEXP paramsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type cors(corsSEXP);
    Rcpp::traits::input_parameter< List >::type init_betas(init_betasSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nj(njSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aj(ajSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(gridCCDr(cors, init_betas, nj, aj, lambdas, params, verbose));
    return __result;
END_RCPP
}
// singleCCDr
List singleCCDr(NumericVector cors, List init_betas, IntegerVector nj, NumericVector aj, double lambda, NumericVector params, int verbose);
RcppExport SEXP ccdr_singleCCDr(SEXP corsSEXP, SEXP init_betasSEXP, SEXP njSEXP, SEXP ajSEXP, SEXP lambdaSEXP, SEXP paramsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type cors(corsSEXP);
    Rcpp::traits::input_parameter< List >::type init_betas(init_betasSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nj(njSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aj(ajSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(singleCCDr(cors, init_betas, nj, aj, lambda, params, verbose));
    return __result;
END_RCPP
}
