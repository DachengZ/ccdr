context("num.edges")

test_that("num.edges works on edgeLists", {
    ### Trivial case
    expect_equal(num.edges(generate_empty_edgeList()), 0)

    ### Non-trivial case
    edgeL <- generate_fixed_edgeList()
    expect_equal(num.edges(edgeL), 5)
})

test_that("num.edges works on SBM", {
    ### Trivial case
    expect_equal(num.edges(generate_empty_SparseBlockMatrixR()), 0)

    ### Non-trivial case
    sbm <- generate_fixed_SparseBlockMatrixR()
    expect_equal(num.edges(sbm), 5)
})

test_that("num.edges works on ccdrFit", {
    ### Trivial case
    cf <- generate_empty_ccdrFit()
    expect_equal(num.edges(cf), 0)

    ### Non-trivial case
    sbm <- generate_fixed_SparseBlockMatrixR()
    cf <- generate_fixed_ccdrFit()
    expect_equal(num.edges(cf), 5)
})

test_that("num.edges works on ccdrPath", {
    ### Trivial case
    cp <- generate_empty_ccdrPath()
    expect_equal(num.edges(cp), rep(0, length(cp)))

    ### Non-trivial case
    cp <- generate_fixed_ccdrPath()
    expect_equal(num.edges(cp), rep(5, length(cp)))
})
