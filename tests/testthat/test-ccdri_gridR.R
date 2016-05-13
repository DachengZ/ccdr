context("ccdri_gridR")

pp <- 10
nn <- 5
X.test <- matrix(runif(pp*pp), ncol = pp)
betas.test <- matrix(runif(pp*pp), ncol = pp)

test_that("ccdri_gridR runs as expected", {
    ### No error / Also tests default values for betas, nlam, gamma, eps, maxIters, alpha
    expect_error(not(ccdri_gridR(X = X.test)))

    ### Also works if X is a data.frame
    expect_error(not(ccdri_gridR(X = as.data.frame(X.test))))
})

test_that("Check input: X", {

    ### Throw error if X has an integer column
    X.test[, 1] <- as.integer(X.test[, 1])
    expect_error(ccdri_gridR(X = X.test))

    ### Check for missing values
    X.test[, 1] <- as.numeric(X.test[, 1]) # re-cast to numeric column
    X.test[1, 1] <- NA
    expect_error(ccdri_gridR(X = X.test))

})

test_that("Check input: lambdas", {

})

test_that("Check input: nlam", {

})

test_that("Check input: rlam", {
    ### eps is non-numeric
    expect_error(ccdri_gridR(X = X.test, rlam = "not a number"))

    ### rlam is negative (note that rlam = 0 is OK)
    expect_error(ccdri_gridR(X = X.test, rlam = -1e-5))

    ### Need to add a test for zero case

})

test_that("Check input: maxIters", {
    ### maxIters is not an integer
    expect_error(ccdri_gridR(X = X.test, maxIters = pi))

    ### maxIters is negative
    expect_error(ccdri_gridR(X = X.test, maxIters = -5))
})

test_that("Check input: alpha", {
    ### alpha is negative
    expect_error(ccdri_gridR(X = X.test, alpha = -5))
})

test_that("Check input: verbose", {
    ### Works when verbose = TRUE
    expect_error(ccdri_gridR(X = X.test, verbose = TRUE))

    ### Works when verbose = FALSE
    expect_error(ccdri_gridR(X = X.test, verbose = FALSE))
})

