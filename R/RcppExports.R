# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

gridCCDr <- function(cors, init_betas, nj, indexj, aj, lambdas, params, verbose) {
    .Call('ccdri_gridCCDr', PACKAGE = 'ccdri', cors, init_betas, nj, indexj, aj, lambdas, params, verbose)
}

singleCCDr <- function(cors, init_betas, nj, indexj, aj, lambda, params, verbose) {
    .Call('ccdri_singleCCDr', PACKAGE = 'ccdri', cors, init_betas, nj, indexj, aj, lambda, params, verbose)
}

