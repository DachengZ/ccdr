#
#  ccdr-main-R.R
#  ccdri
#
#  Created by Bryon Aragam (local) on 2/4/15.
#  Copyright (c) 2014-2015 Bryon Aragam (local). All rights reserved.
#

#
# PACKAGE CCDR: Main CCDr methods
#
#   CONTENTS:
#     ccdr.run
#     ccdri_call
#     ccdri_gridR
#     ccdri_singleR
#

###--- These two lines are necessary to import the auto-generated Rcpp methods in RcppExports.R---###
#' @useDynLib ccdri
#' @importFrom Rcpp sourceCpp
NULL

#' Main CCDr Algorithm
#'
#' Estimate a Bayesian network (directed acyclic graph) from observational data using the
#' CCDr algorithm as described in \href{http://arxiv.org/abs/1401.0852}{Aragam and Zhou (2015), JMLR}.
#'
#' Instead of producing a single estimate, this algorithm computes a solution path of estimates based
#' on the values supplied to \code{lambdas} or \code{lambdas.length}. The CCDr algorithm approximates
#' the solution to a nonconvex optimization problem using coordinate descent. Instead of AIC or BIC,
#' CCDr uses continuous regularization based on concave penalties such as the minimax concave penalty
#' (MCP).
#'
#' This implementation includes two options for the penalty: (1) MCP, and (2) L1 (or Lasso). This option
#' is controlled by the \code{gamma} argument.
#'
#' @param data Data matrix. Must be numeric and contain no missing values.
#' @param betas Initial guess for the algorithm. Represents the weighted adjacency matrix
#'              of a DAG where the algorithm will begin searching for an optimal structure.
#' @param lambdas (optional) Numeric vector containing a grid of lambda values (i.e. regularization
#'                parameters) to use in the solution path. If missing, a default grid of values will be
#'                used based on a decreasing log-scale  (see also \link{generate.lambdas}).
#' @param lambdas.length Integer number of values to include in the solution path. If \code{lambdas}
#'                       has also been specified, this value will be ignored. Note also that the final
#'                       solution path may contain fewer estimates (see
#'                       \code{alpha}).
#' @param gamma Value of concavity parameter. If \code{gamma > 0}, then the MCP will be used
#'              with \code{gamma} as the concavity parameter. If \code{gamma < 0}, then the L1 penalty
#'              will be used and this value is otherwise ignored.
#' @param error.tol Error tolerance for the algorithm, used to test for convergence.
#' @param max.iters Maximum number of iterations for each internal sweep.
#' @param alpha Threshold parameter used to terminate the algorithm whenever the number of edges in the
#'              current estimation is \code{> alpha * ncol(data)}.
#' @param verbose \code{TRUE / FALSE} whether or not to print out progress and summary reports.
#'
#' @return A \code{\link{ccdrPath-class}} object.
#'
#' @examples
#'
#' \dontrun{
#'
#' ### Generate some random data
#' dat <- matrix(rnorm(1000), nrow = 20)
#'
#' ### Run with default settings
#' ccdr.run(data = dat)
#'
#' ### Optional: Adjust settings
#' pp <- ncol(dat)
#' init.betas <- matrix(0, nrow = pp, ncol = pp)              # initialize algorithm with a random initial value
#' init.betas[1,2] <- init.betas[1,3] <- init.betas[4,2] <- 1 #
#' ccdr.run(data = dat, betas = init.betas, lambdas.length = 10, alpha = 10, verbose = TRUE)
#' }
#'
#' @export
ccdr.run <- function(data,
                     intervention = NULL, ## now use a seperate vector indicating which node (1, ..., pp) is fixed; pp+1 means no intervention
                                             ## sort? maybe add another wrapper
                     betas,
                     lambdas,
                     lambdas.length = NULL,
                     gamma = 2.0,
                     error.tol = 1e-4,
                     max.iters = NULL,
                     alpha = 10,
                     verbose = FALSE
){
    ### This is just a wrapper for the internal implementation given by ccdri_call
    ccdri_call(data = data,
              intervention = intervention,
              betas = betas,
              lambdas = lambdas,
              lambdas.length = lambdas.length,
              gamma = gamma,
              error.tol = error.tol,
              rlam = NULL,
              max.iters = max.iters,
              alpha = alpha,
              verbose = verbose)
} # END CCDR.RUN

# ccdri_call
#
#   Handles most of the bookkeeping for CCDr. Sets default values and prepares arguments for
#    passing to ccdri_gridR and ccdri_singleR. Some type-checking as well, although most of
#    this is handled internally by ccdri_gridR and ccdri_singleR.
#
ccdri_call <- function(data,
                      intervention,
                      betas,
                      lambdas,
                      lambdas.length,
                      gamma,
                      error.tol,
                      rlam,
                      max.iters,
                      alpha,
                      verbose = FALSE
){
    ### Check data
    if(!check_if_data_matrix(data)) stop("Data must be either a data.frame or a numeric matrix!")
    if(count_nas(data) > 0) stop(paste0(count_nas(data), " missing values detected!"))

    ### Get the dimensions of the data matrix
    nn <- nrow(data)
    pp <- ncol(data)

    ## check intervention
    if(!is.null(intervention)) {
        if(length(intervention) != nn) stop("Intervention size does not match.")
        if(any(!is.integer(intervention))) stop("Intervention labels are not integers.")
        if(any(intervention <= 0) || any(intervention > pp + 1)) stop("Intervention out of range.")
    } else intervention <- rep(pp + 1, nn)

    nj <- rep(0, pp)
    for(j in 1:pp) { ## include 0 here or not?
        nj[j] <- sum(intervention != j) ## optimize for sorted column?
    }

    ### Use default values for lambda if not specified
    if(missing(lambdas)){
        if(is.null(lambdas.length)){
            stop("Both lambdas and lambdas.length unspecified: Must specify a value for at least one of these arguments!")
        } else{
            ### Check lambdas.length if specified
            if(!is.numeric(lambdas.length)) stop("lambdas.length must be numeric!")
            if(lambdas.length <= 0) stop("lambdas.length must be positive!")
        }

        if(missing(rlam)){
            ### Even though ccdri_call should never be called on its own, this behaviour is left for testing backwards-compatibility
            stop("rlam must be specified if lambdas is not explicitly specified.")
        } else if(is.null(rlam)){
            ### rlam = NULL is used as a sentinel value to indicate a default value should be used
            rlam <- 1e-2
        } else{
            ### Check rlam if specified
            if(!is.numeric(rlam)) stop("rlam must be numeric!")
            if(rlam < 0) stop("rlam must be >= 0!")
        }

        # If no grid of lambdas is passed, then use the standard log-scale that starts at
        #  max.lam = sqrt(nn) and descends to min.lam = rlam * max.lam
        lambdas <- generate.lambdas(lambda.max = sqrt(nn),
                                    lambdas.ratio = rlam,
                                    lambdas.length = as.integer(lambdas.length),
                                    scale = "log")
    }

    ### Check lambdas
    if(!is.numeric(lambdas)) stop("lambdas must be a numeric vector!")
    if(any(lambdas < 0)) stop("lambdas must contain only nonnegative values!")

#     if(length(lambdas) != lambdas.length){
#         warning("Length of lambdas vector does not match lambdas.length. The specified lambdas vector will be used and lambdas.length will be overwritten.")
#     }

    ### By default, set the initial guess for betas to be all zeroes
    if(missing(betas)){
        betas <- matrix(0, nrow = pp, ncol = pp)
        # betas <- SparseBlockMatrixR(betas) # 2015-03-26: Deprecated and replaced with .init_sbm below
        betas <- .init_sbm(betas, rep(0, pp))

        # If the initial matrix is the zero matrix, indexing does not matter so we don't need to use reIndexC here
        #   Still need to set start = 0, though.
        betas$start <- 0
    } # Type-checking for betas happens in ccdri_singleR

    # This parameter can be set by the user, but in order to prevent the algorithm from taking too long to run
    #  it is a good idea to keep the threshold used by default which is O(sqrt(pp))
    if(is.null(max.iters)){
        max.iters <- 2 * max(10, sqrt(pp))
    }

    t1.cor <- proc.time()[3]
    #     cors <- cor(data)
    #     cors <- cors[upper.tri(cors, diag = TRUE)]
    corlist <- cor_vector_intervention(data, intervention) ## cors <- cor_vector(data)
    cors <- corlist$cors
    indexj <- corlist$indexj
    t2.cor <- proc.time()[3]

    fit <- ccdri_gridR(cors = cors, ## now a longer vector
                      # data = data,
                      as.integer(pp),
                      as.integer(nn), ## can we replace nn with nj now?
                      as.integer(nj), ## added nj
                      as.integer(indexj),
                      betas,
                      as.numeric(lambdas),
                      as.numeric(gamma),
                      as.numeric(error.tol),
                      as.integer(max.iters),
                      as.numeric(alpha),
                      verbose)

    fit <- lapply(fit, ccdrFit.list)    # convert everything to ccdrFit objects
    ccdrPath.list(fit)                  # wrap as ccdrPath object
} # END ccdri_CALL

# ccdri_gridR
#
#   Main subroutine for running the CCDr algorithm on a grid of lambda values.
ccdri_gridR <- function(cors, ## now a longer vector
                       # data,
                       pp, nn, ##
                       nj, ## added nj
                       indexj,
                       betas,
                       lambdas,
                       gamma,
                       eps,
                       maxIters,
                       alpha,
                       verbose
){

    ### Check alpha
    if(!is.numeric(alpha)) stop("alpha must be numeric!")
    if(alpha < 0) stop("alpha must be >= 0!")

    ### nlam is now set automatically
    nlam <- length(lambdas)

    ccdr.out <- list()
    for(i in 1:nlam){
        if(verbose) message("Working on lambda = ", round(lambdas[i], 5), " [", i, "/", nlam, "]")

        t1.ccdr <- proc.time()[3]
        ccdr.out[[i]] <- ccdri_singleR(cors = cors, ## now a longer vector
                                      pp, nn, ##
                                      nj, ## now added nj
                                      indexj,
                                      betas,
                                      lambdas[i],
                                      gamma = gamma,
                                      eps = eps,
                                      maxIters = maxIters,
                                      alpha = alpha,
                                      verbose = verbose
        )
        t2.ccdr <- proc.time()[3]

        sbm <- ccdr.out[[i]]$sbm
        betas <- reIndexC(sbm) # use C-friendly indexing

        if(verbose){
            test.nedge <- sum(as.matrix(betas) != 0)
            message("  Estimated number of edges: ", ccdr.out[[i]]$nedge, " / test = ", test.nedge)
            # message("  Estimated total variance: ", sum(1 / (betas$sigmas)^2))
        }


        # # ## calculate log-likelihood
        # sbm <- to_B(sbm)
        # B <- as.matrix(sbm)
        # s <- sbm$sigmas
        #
        # intervention <- data[, pp + 1]
        # logl <- 0
        # for(j in 1:pp) {
        #     loglj <- sum((data[intervention != j, j] - data[intervention != j, -c(j, pp + 1), drop = FALSE] %*% B[-j, j, drop = FALSE])^2)
        #     loglj <- loglj / 2 / s[j] + nj[j] * log(s[j]) / 2
        #     logl <- logl + loglj
        # }
        # ccdr.out[[i]]$logl <- logl

        # 7-16-14: Added code below to check edge threshold via alpha parameter
        if(ccdr.out[[i]]$nedge > alpha * pp){
            if(verbose) message("Edge threshold met, terminating algorithm with ", ccdr.out[[i-1]]$nedge, " edges.")
            break
        }
    }

    ccdr.out[1:(i-1)] # only return up to i - 1 since the last (ith) model would not have finished running anyway
} # END ccdri_GRIDR

# ccdri_singleR
#
#   Internal subroutine for handling calls to singleCCDr: This is the only place where C++ is directly
#    called. Type-checking is strongly enforced here.
ccdri_singleR <- function(cors, ## now a longer vector
                         pp, nn, ## remove nn?
                         nj, ## added nj
                         indexj,
                         betas,
                         lambda,
                         gamma,
                         eps,
                         maxIters,
                         alpha,     # 2-9-15: No longer necessary in ccdri_singleR, but needed since the C++ call asks for it
                         verbose = FALSE
){

    ### Check cors
    if(!is.numeric(cors)) stop("cors must be a numeric vector!")
    if(length(cors) != length(unique(indexj))*pp*(pp+1)/2) stop(paste0("cors has incorrect length: Expected length = ", pp*pp*(pp+1)/2, " input length = ", length(cors)))
    ## longer due to intervention

    ### Check dimension parameters
    if(!is.integer(pp) || !is.integer(nn)) stop("Both pp and nnmust be integers!")
    if(pp <= 0 || nn <= 0) stop("Both pp and nn must be positive!")

    ## check nj
    if(any(!is.integer(nj)) || any(nj < 0)) stop("nj must be non-negative integers!")
    if(length(nj) != pp) stop(paste0("nj size does not match:", paste(nj, collapse = ", "), ". Check it out."))
    ## check sum(nj) - (num of non-intervention) == (pp-1) * nn
    aj <- nj / nn

    ### Check betas
    if(check_if_matrix(betas)){ # if the input is a matrix, convert to SBM object
        betas <- SparseBlockMatrixR(betas) # if betas is non-numeric, SparseBlockMatrixR constructor should throw error
        betas <- reIndexC(betas) # use C-friendly indexing
    } else if(!is.SparseBlockMatrixR(betas)){ # otherwise check that it is an object of class SparseBlockMatrixR
        stop("Incompatible data passed for betas parameter: Should be either matrix or list in SparseBlockMatrixR format.")
    }

    ### Check lambda
    if(!is.numeric(lambda)) stop("lambda must be numeric!")
    if(lambda < 0) stop("lambda must be >= 0!")

    ### Check gamma
    if(!is.numeric(gamma)) stop("gamma must be numeric!")
    if(gamma < 0 && gamma != -1) stop("gamma must be >= 0 (MCP) or = -1 (Lasso)!")

    ### Check eps
    if(!is.numeric(eps)) stop("eps must be numeric!")
    if(eps <= 0){
        if(eps < 0) stop("eps must be >= 0!")
        if(eps == 0) warning("eps is set to zero: This may cause the algorithm to fail to converge, and maxIters will be used to terminate the algorithm.")
    }

    ### Check maxIters
    if(!is.integer(maxIters)) stop("maxIters must be an integer!")
    if(maxIters <= 0) stop("maxIters must be > 0!")

    ### alpha check is in ccdri_gridR

    # if(verbose) cat("Opening C++ connection...")
    t1.ccdr <- proc.time()[3]
    ccdr.out <- singleCCDr(cors,
                           betas,
                           nj, ## now added nj
                           indexj,
                           aj,
                           lambda,
                           c(gamma, eps, maxIters, alpha),
                           verbose = verbose)
    t2.ccdr <- proc.time()[3]
    # if(verbose) cat("C++ connection closed. Total time in C++: ", t2.ccdr-t1.ccdr, "\n")

    #
    # Convert output back to SBM format
    #
    ccdr.out <- list(sbm = SparseBlockMatrixR(list(rows = ccdr.out$rows, vals = ccdr.out$vals, blocks = ccdr.out$blocks, sigmas = ccdr.out$sigmas, start = 0)),
                     lambda = ccdr.out$lambda,
                     nedge = ccdr.out$length,
                     pp = pp,
                     nn = nn, ## do we still need to include nj here?
                     time = t2.ccdr - t1.ccdr
                     )
    ccdr.out$sbm <- reIndexR(ccdr.out$sbm)

    # ccdrFit(ccdr.out)
    ccdr.out
} # END ccdri_SINGLER
