//
//  rcpp_wrap.cpp
//  ccdr_proj
//
//  Created by Bryon Aragam on 3/24/14.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

#ifndef rcpp_wrap_h
#define rcpp_wrap_h

#include "defines.h"

#ifdef _COMPILE_FOR_RCPP_

#include <Rcpp.h>
#include "algorithm.h"

using namespace Rcpp;

//------------------------------------------------------------------------------/
//   WRAPPER FUNCTIONS FOR R / RCPP INTEGRATION
//------------------------------------------------------------------------------/

//
// These functions provide the necessary wrapper functions for providing compatibility with Rcpp and R.
//   By keeping these functions in a separate header, we can ensure that the Rcpp dependencies are only
//   loaded when compiling from R---the other files in the codebase should be completely independent of
//   Rcpp.
//
// One exception to this is algorithm.h, but these dependencies are presently safely handled by using the
//   _COMPILE_FOR_RCPP_ ifdef.
//

//
// Translation key for C++ to Rcpp
//
//   Variable Types
//      std::vector<double>: NumericVector
//      std::vector<int>: IntegerVector
//
//   Functions
//      as<>: convert Rcpp object to C++ object
//      wrap<>: convert C++ to Rcpp object (type handled automatically)
//

// [[Rcpp::export]]
List gridCCDr(NumericVector cors,
              List init_betas,
              unsigned int nn,
              IntegerVector nj, // added nj
              NumericVector lambdas,
              NumericVector params,
              int verbose
              ){
    SparseBlockMatrix betas = SparseBlockMatrix(init_betas);

    #ifdef _DEBUG_ON_
        //
        // log.h logging
        //
        FILE* pFile = fopen("/Users/Dacheng/Desktop/ccdr_proj_LOG_FILE.txt", "w");
        Output2FILE::Stream() = pFile;
        FILELog::ReportingLevel() = logDEBUG1;

        FILE_LOG(logINFO) << "Log file opened.";
    #endif

    std::vector<SparseBlockMatrix> grid_betas;
    grid_betas = gridCCDr(as< std::vector<double> >(cors),
                          betas,
                          nn,
                          as< std::vector<int> >(nj), // added nj
                          as< std::vector<double> >(lambdas),
                          as< std::vector<double> >(params),
                          verbose);

    std::vector<List> return_betas;
    for(int i = 0; i < grid_betas.size(); ++i){
        return_betas.push_back(grid_betas[i].get_R(lambdas[i]));
    }

    return wrap(return_betas);
}

// [[Rcpp::export]]
List singleCCDr(NumericVector cors,
                List init_betas,
                IntegerVector nj, // added nj
                double lambda,
                NumericVector params,
                int verbose
                ){

    SparseBlockMatrix betas = SparseBlockMatrix(init_betas);

    betas = singleCCDr(as< std::vector<double> >(cors),
                       betas,
                       //nn,
                       as< std::vector<int> >(nj), // added nj
                       lambda,
                       as< std::vector<double> >(params),
                       verbose);
    //
    // Need to manually recompute active set size when calling singleCCDr directly from R,
    //   as opposed to within gridCCDr, which automatically recomputes the active set size
    //
    betas.recomputeActiveSetSize(true);


    return betas.get_R(lambda);
}

//---------------------------------------------------------------------------------------------------//
// ***IF THIS CODE THROWS ANY ERRORS, MOVE THIS DEFINITION BACK TO THE END OF SparseBlockMatrix.h***
//
//
// Returns the SparseBlockMatrix as an R list (using Rcpp); for passing data back to R
//  Two cases:
//  1) Include lambda in list (lambda_R >= 0)
//  2) Ignore lambda (lambda_R < 0)
List SparseBlockMatrix::get_R(double lambda_R){
    if(lambda_R < 0)
        return List::create(_["rows"] = wrap(rows), _["vals"] = wrap(vals), _["sigmas"] = wrap(sigmas), _["blocks"] = wrap(blocks), _["length"] = wrap(activeSetLength));
    else
        return List::create(_["rows"] = wrap(rows), _["vals"] = wrap(vals), _["sigmas"] = wrap(sigmas), _["blocks"] = wrap(blocks), _["length"] = wrap(activeSetLength), _["lambda"] = wrap(lambda_R));
}
//---------------------------------------------------------------------------------------------------//

#endif

#endif
