#
#  s3-SparseBlockMatrixR.R
#  ccdri
#
#  Created by Bryon Aragam (local) on 2/4/15.
#  Copyright (c) 2014-2015 Bryon Aragam (local). All rights reserved.
#

#------------------------------------------------------------------------------#
# SparseBlockMatrixR S3 Class for R
#------------------------------------------------------------------------------#

#
# SparseBlockMatrixR S3 class skeleton
#
# Data
# * list rows
# * list vals
# * list blocks
# * numeric sigmas
# * integer start
#
# Methods
# * is.SparseBlockMatrixR
# * reIndexC.SparseBlockMatrixR
# * reIndexR.SparseBlockMatrixR
# * SparseBlockMatrixR.list
# * SparseBlockMatrixR.sparse
# * SparseBlockMatrixR.matrix
# * as.SparseBlockMatrixR.list
# * as.SparseBlockMatrixR.sparse
# * as.SparseBlockMatrixR.matrix
# * as.list.SparseBlockMatrixR
# * as.matrix.SparseBlockMatrixR
# * as.edgeList.SparseBlockMatrixR
# * get.adjacency.matrix.SparseBlockMatrixR
# * num.nodes.SparseBlockMatrixR
# * num.edges.SparseBlockMatrixR
# * is.zero.SparseBlockMatrixR
# * .init_sbm
# * to_B.SparseBlockMatrixR
#

#
# A convenience class to make easier sharing data between R and C++ easier. This class mimics the structure
#   of the C++ class 'SparseBlockMatrix' (note the name difference to differentiate the two) as a list in R,
#   which makes it easy to use Rcpp to pass a sparse structure between R and C++. This class is NOT intended
#   to be a general purpose sparse data structure for R; instead, it simply streamlines the connection between
#   R and C++.
#
# A SparseBlockMatrixR object consists of four main components:
#   1) rows - a list of integer vectors
#   2) vals - a list of numeric vectors
#   3) blocks - a list of integer vectors
#   4) sigmas - a numeric vector
#
# There is also a fifth component which identifies whether the indexing in the object begins at 0 or 1.
#   This is needed for bookkeeping and ensuring coherent translation between R and C++.
#   5) start - 0 or 1
#
# These components all exactly reflect their purpose in SparseBlockMatrix.h; see that file for more details.
#
#   NOTES:
#       1) Since C++ begins indexing at 0, we have included two functions for re-indexing between the two conventions,
#          defined in reIndexR and reIndexC. The 'start' flag ensures we keep track of where indexing begins.
#
#

#------------------------------------------------------------------------------#
# is.SparseBlockMatrixR
#
#' @export
is.SparseBlockMatrixR <- function(sbm){
    inherits(sbm, "SparseBlockMatrixR")
} # END IS.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# reIndexC.SparseBlockMatrixR
#  Re-indexing TO C for SparseBlockMatrixR objects
#
reIndexC.SparseBlockMatrixR <- function(sbm){
    #
    # Using lapply does NOT work here: if one of the list elements is an empty vector, adding 1 will
    #  mysteriously coerce it to a numeric vector (should be integer!). Not sure why this happens, but
    #  we shouldn't allow R to secretly change our data types without our permission. Use a for loop
    #  instead.
    #
    # UPDATE 05/13/14: Using '1L' keeps everything as an integer. This makes sense since '1' is a numeric
    #                   literal in R, while the corresponding literal for an integer is '1L'.
    #
    sbm$rows <- lapply(sbm$rows, function(x){ x - 1L})
    if(length(sbm$blocks) > 0) sbm$blocks <- lapply(sbm$blocks, function(x){ x - 1L})

    sbm$start <- 0

    sbm
} # END REINDEXC.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# reIndexR.SparseBlockMatrixR
#  Re-indexing TO R for SparseBlockMatrixR objects
#
reIndexR.SparseBlockMatrixR <- function(sbm){
    #
    # Using lapply does NOT work here: if one of the list elements is an empty vector, adding 1 will
    #  coerce it to a numeric vector (should be integer!). Using '1L' keeps everything as an integer.
    #  This makes sense since '1' is a numeric literal in R, while the corresponding literal for an
    #  integer is '1L'.
    #
    sbm$rows <- lapply(sbm$rows, function(x){ x + 1L})
    if(length(sbm$blocks) > 0) sbm$blocks <- lapply(sbm$blocks, function(x){ x + 1L})

    sbm$start <- 1

    sbm
} # END REINDEXR.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.list
#  List constructor
#
SparseBlockMatrixR.list <- function(li){

    if( !is.list(li)){
        stop("Input must be a list!")
    }

    if( length(li) != 5 || !setequal(names(li), c("rows", "vals", "blocks", "sigmas", "start"))){
        stop("Input is not coercable to an object of type SparseBlockMatrixR, check list for the following elements: rows (list), vals (list), blocks (list), sigmas (numeric), start (integer)")
    }

    if(!is.list(li$rows) || !is.list(li$vals)){
        stop("rows and vals must both be lists of length pp!")
    }

    if( length(li$rows) != length(li$vals)){
        #
        # We enforce that rows and vals must have the same length, but relax this assumption for blocks and sigmas
        #  since they are mostly internal to the CCDr algorithm, and once the algorithm has run we may want to free
        #  up the memory associated with these elements
        #
        stop("rows and vals have different sizes; should all have the same length (pp)!!")
    }

    structure(li, class = "SparseBlockMatrixR")
} # END SPARSEBLOCKMATRIXR.LIST

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.sparse
#  sparse object constructor
#
SparseBlockMatrixR.sparse <- function(sp){

    if( !is.sparse(sp)){
        stop("Input must be a sparse object!")
    } else if(sp$dim[1] != sp$dim[2]){
        stop("Input must be square!")
    }

    pp <- sp$dim[1]
    if(sp$start == 0) sp <- reIndexR(sp) # re-index rows and cols to start at 1 if necessary

    sbm.rows <- vector("list", length = pp)
    sbm.vals <- vector("list", length = pp)
    sbm.blocks <- vector("list", length = pp)

    warning("Attempting to coerce sparse object to SparseBlockMatrixR with no data for sigmas: \n   Setting sigma_j = 0 by default.")
    sbm.sigmas <- rep(0, pp) ### 2015-03-25: check this default!

    # how to vectorize this???
    for(j in 1:pp){
        # Clear out the jth entry in the lists to be an empty vector
        sbm.rows[[j]] <- integer(0)
        sbm.vals[[j]] <- numeric(0)
        sbm.blocks[[j]] <- integer(0)
    }

    for(j in 1:pp){

        thisColIdx <- which(sp$cols == j)
        rows <- as.integer(sp$rows[thisColIdx])
        for(k in seq_along(rows)){
            row <- rows[k]

            sbm.rows[[j]] <- c(sbm.rows[[j]], row)
            sbm.rows[[row]] <- c(sbm.rows[[row]], j)

            sbm.vals[[j]] <- c(sbm.vals[[j]], sp$vals[thisColIdx[k]])
            sbm.vals[[row]] <- c(sbm.vals[[row]], 0)

            # vals[rows[j][k]][block[j][k]] = beta_ji
            sbm.blocks[[j]] <- c(sbm.blocks[[j]], length(sbm.rows[[row]]))
            sbm.blocks[[row]] <- c(sbm.blocks[[row]], length(sbm.rows[[j]]))
        }

    }

    names(sbm.rows) <- names(sbm.vals) <- names(sbm.blocks) <- as.character(1:pp)

    #
    # NOTE: We use R-indexing by default. This can be changed by using reIndexC if necessary.
    #
    SparseBlockMatrixR.list(list(rows = sbm.rows, vals = sbm.vals, blocks = sbm.blocks, sigmas = sbm.sigmas, start = 1))
} # END SPARSEBLOCKMATRIXR.SPARSE

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.matrix
#  matrix constructor
#
SparseBlockMatrixR.matrix <- function(m){

    if(nrow(m) != ncol(m)) stop("Input matrix must be square!")

    sp <- as.sparse(m)

    SparseBlockMatrixR.sparse(sp)
} # END SPARSEBLOCKMATRIXR.MATRIX

#------------------------------------------------------------------------------#
# as.SparseBlockMatrixR.list
#  Convert FROM list TO SparseBlockMatrixR
#
as.SparseBlockMatrixR.list <- function(li){
    SparseBlockMatrixR.list(li)
} # END AS.SPARSEBLOCKMATRIXR.LIST

#------------------------------------------------------------------------------#
# as.SparseBlockMatrixR.sparse
#  Convert FROM sparse TO SparseBlockMatrixR
#
#' @export
as.SparseBlockMatrixR.sparse <- function(sp){
    SparseBlockMatrixR.sparse(sp)
} # END AS.SPARSEBLOCKMATRIXR.SPARSE

#------------------------------------------------------------------------------#
# as.SparseBlockMatrixR.matrix
#  Convert FROM matrix TO SparseBlockMatrixR
#
#' @export
as.SparseBlockMatrixR.matrix <- function(m){
    SparseBlockMatrixR.matrix(m)
} # END AS.SPARSEBLOCKMATRIXR.MATRIX

#------------------------------------------------------------------------------#
# as.list.SparseBlockMatrixR
#  Convert FROM SparseBlockMatrixR TO list
#  Even though internally the SBM object is a list, we must still manually define this function
#
#' @export
as.list.SparseBlockMatrixR <- function(sbm){
    list(rows = sbm$rows, vals = sbm$vals, blocks = sbm$blocks, sigmas = sbm$sigmas, start = sbm$start)
} # END AS.LIST.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# as.matrix.SparseBlockMatrixR
#  Convert FROM SparseBlockMatrixR TO matrix
#
#' @export
as.matrix.SparseBlockMatrixR <- function(sbm){
    pp <- length(sbm$rows)
    m <- matrix(0, nrow = pp, ncol = pp)

    ### 2015-03-02: Why was I using diag to construct this matrix?
    # m <- diag(rep(0, pp))

    if(sbm$start == 0) sbm <- reIndexR(sbm)

    for(j in 1:pp){
        m[sbm$rows[[j]], j] <- sbm$vals[[j]]
    }

    ### 2015-03-02: Do not need to set dim attribute of matrix! (Already set by default constructor)
    # attributes(m)$dim <- c(pp, pp)
    # attributes(m)$dimnames <- list()
    rownames(m) <- as.character(1:nrow(m))
    colnames(m) <- as.character(1:ncol(m))

    m
} # END AS.MATRIX.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# as.edgeList.SparseBlockMatrixR
# Coerce SBM to edge list
#
#' @export
as.edgeList.SparseBlockMatrixR <- function(sbm){
    #
    # We have to be careful in obtaining the edge list of a SparseBlockMatrixR object:
    #  It is NOT the same as the rows slot since some of these components may have
    #  zero edge weights (see docs for SparseBlockMatrixR for explanation). Thus, in
    #  order to obtain the edge list, we need to check which indices in the rows slot
    #  have nonzero edge weights.
    #
    # y = rows, x = vals : Select the elements of rows which have nonzero values in vals,
    #                       accouting for possible round-off (hence .MACHINE_EPS).
    #
    el <- mapply(function(x, y){ y[which(abs(x) > .MACHINE_EPS)]}, sbm$vals, sbm$rows)

    edgeList.list(el)
} # AS.EDGELIST.SPARSEBLOCKMATRIXR

#' @export
#' @describeIn get.adjacency.matrix Convert internal \code{SparseBlockMatrixR} representation to an adjacency matrix
get.adjacency.matrix.SparseBlockMatrixR <- function(sbm){
    get.adjacency.matrix.edgeList(as.edgeList.SparseBlockMatrixR(sbm))
} # END GET.ADJACENCY.MATRIX.SPARSEBLOCKMATRIXR


#' @export
#' @describeIn num.nodes
num.nodes.SparseBlockMatrixR <- function(sbm){
    ### The number of nodes should be exactly the same as the length of the rows list
    length(sbm$rows)
} # END NUM.NODES.SPARSEBLOCKMATRIXR

#' @export
#' @describeIn num.edges
num.edges.SparseBlockMatrixR <- function(sbm){
    ### The number of nodes should be exactly the same as the length of the rows list
    num.edges(as.edgeList.SparseBlockMatrixR(sbm))
} # END NUM.EDGES.SPARSEBLOCKMATRIXR

#' @export
#' @describeIn is.zero
is.zero.SparseBlockMatrixR <- function(x){
    check_if_zero <- (length(unlist(x$sbm$rows)) == 0)

    check_if_zero
} # END IS.ZERO.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# .init_sbm
# Internal function for initializing a SparseBlockMatrixR object directly
#  from a matrix AND a sigmas vector
#
.init_sbm <- function(init_matrix, init_sigmas){
    stopifnot(check_if_matrix(init_matrix))
    stopifnot(nrow(init_matrix) == ncol(init_matrix))

    stopifnot(is.numeric(init_sigmas))
    stopifnot(length(init_sigmas) == nrow(init_matrix))

    sbm <- suppressWarnings(SparseBlockMatrixR(init_matrix)) # suppress warnings since we are working in a controlled environment
    sbm$sigmas <- init_sigmas

    sbm
} # END .INIT_SBM

#------------------------------------------------------------------------------#
# .to_B.SparseBlockMatrixR
# Internal function to convert estimates from the (Rho, R) parametrization to
#  the standard (B, Omega) parametrization.
#
to_B.SparseBlockMatrixR <- function(sbm){
    ### Need to re-parametrize the betas FIRST
    sbm$vals <- mapply(function(x, y) x/y, sbm$vals, sbm$sigmas) # Divide each vals vector by sigmas[j]
    sbm$sigmas <- 1/ (sbm$sigmas)^2

    sbm
}
