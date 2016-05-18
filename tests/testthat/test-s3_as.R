context("sparse / SparseBlockMatrixR to matrix conversion")

#
# NOTE: Need to use is_equivalent_to in order to ignore attributes!
#

### matrix -> sparse
test_that("as.matrix -> as.sparse -> as.matrix makes no changes", {

    ### Test the zero matrix
    m <- matrix(rep(0, 1), ncol = 1)
    expect_that(as.matrix(as.sparse(as.matrix(m))), is_equivalent_to(m))
    m <- matrix(rep(0, 4), ncol = 2)
    expect_that(as.matrix(as.sparse(as.matrix(m))), is_equivalent_to(m))

    ### 2-27-2015: sparse currently not allowed to handle non-square matrices: Why though?
#     m <- matrix(rep(0, 12), ncol = 3) # non-square case
#     expect_that(as.matrix(as.sparse(as.matrix(m))), is_equivalent_to(m))

    ### Test a randomly generated matrix
    m <- random.sparse(5, 5)
    expect_that(as.matrix(as.sparse(as.matrix(m))), is_equivalent_to(m))

    ### Test a randomly generated DAG
    m <- random.dag.matrix(10, 10)
    expect_that(as.matrix(as.sparse(as.matrix(m))), is_equivalent_to(m))
})

### matrix -> SparseBlockMatrixR
test_that("as.matrix -> as.SparseBlockMatrixR -> as.matrix makes no changes", {

    ### Test the zero matrix
    m <- matrix(rep(0, 1), ncol = 1)
    expect_that(as.matrix(as.SparseBlockMatrixR(as.matrix(m))), is_equivalent_to(m))
    m <- matrix(rep(0, 4), ncol = 2)
    expect_that(as.matrix(as.SparseBlockMatrixR(as.matrix(m))), is_equivalent_to(m))

    ### NOTE: Cannot test on random sparse matrix since SBM class ASSUMES a block structure,
    ###        i.e. induced by a DAG

    ### Test a randomly generated DAG
    m <- random.dag.matrix(10, 10)
    expect_that(as.matrix(as.SparseBlockMatrixR(as.matrix(m))), is_equivalent_to(m))
})

### sparse -> SparseBlockMatrixR
test_that("as.sparse -> as.SparseBlockMatrixR -> as.sparse makes no changes", {

    ### Test the zero matrix
    m <- matrix(rep(0, 1), ncol = 1)
    sp <- as.sparse(m)
    expect_that(as.sparse(as.SparseBlockMatrixR(sp)), equals(sp))

    m <- matrix(rep(0, 4), ncol = 2)
    sp <- as.sparse(m)
    expect_that(as.sparse(as.SparseBlockMatrixR(sp)), equals(sp))

    ### NOTE: Cannot test on random sparse matrix since SBM class ASSUMES a block structure,
    ###        i.e. induced by a DAG

    ### Test a randomly generated DAG
    m <- random.dag.matrix(10, 10)
    sp <- as.sparse(m)
    expect_that(as.sparse(as.SparseBlockMatrixR(sp)), equals(sp))
})

### SparseBlockMatrixR -> edgeList -> matrix
test_that("as.SparseBlockMatrixR -> as.edgeList -> as.matrix makes no changes", {

    ### Test the zero matrix
    m <- matrix(rep(0, 1), ncol = 1)
    sbm <- as.SparseBlockMatrixR(m)
    m.sbm <- as.matrix(sbm)
    m.edgeL <- as.matrix(as.edgeList.SparseBlockMatrixR(sbm))
    expect_that(m.edgeL, is_equivalent_to(m.sbm))

    m <- matrix(rep(0, 4), ncol = 2)
    sbm <- as.SparseBlockMatrixR(m)
    m.sbm <- as.matrix(sbm)
    m.edgeL <- as.matrix(as.edgeList.SparseBlockMatrixR(sbm))
    expect_that(m.edgeL, is_equivalent_to(m.sbm))

    ### NOTE: Cannot test on random sparse matrix since SBM class ASSUMES a block structure,
    ###        i.e. induced by a DAG

    ### Test a randomly generated DAG
    m <- random.dag.matrix(10, 10)
    sbm <- as.SparseBlockMatrixR(m)
    m.sbm <- as.matrix(sbm)
    m.edgeL <- as.matrix(as.edgeList.SparseBlockMatrixR(sbm))

    m.sbm[m.sbm != 0] <- 1 # coerce to unweighted adjacencies
    expect_that(m.edgeL, is_equivalent_to(m.sbm))
})
