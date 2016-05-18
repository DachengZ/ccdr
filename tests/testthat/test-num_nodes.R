context("num.nodes")

test_that("num.nodes works on edgeLists", {
    ### Trivial case
    expect_equal(num.nodes(generate_empty_edgeList()), 0)

    ### Non-trivial case
    edgeL <- generate_fixed_edgeList()
    expect_equal(num.nodes(edgeL), 5)
})

test_that("num.nodes works on SBM", {
    ### Trivial case
    expect_equal(num.nodes(generate_empty_SparseBlockMatrixR()), 0)

    ### Non-trivial case
    sbm <- generate_fixed_SparseBlockMatrixR()
    expect_equal(num.nodes(sbm), 5)
})

test_that("num.nodes works on ccdrFit", {
    ### Trivial case
    cf <- generate_empty_ccdrFit()
    expect_equal(num.nodes(cf), 0)

    ### Non-trivial case
    cf <- generate_fixed_ccdrFit()
    expect_equal(num.nodes(cf), 5)
})

test_that("num.nodes works on ccdrPath", {
    ### Trivial case
    cp <- generate_empty_ccdrPath()
    expect_equal(num.nodes(cp), 0)

    ### Non-trivial case
    cp <- generate_fixed_ccdrPath()
    expect_equal(num.nodes(cp), 5)
})
