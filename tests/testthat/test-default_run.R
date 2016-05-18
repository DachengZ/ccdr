context("Default run")

library("pcalg")

### Set test parameters
pp <- 10
nn <- 100
ss <- 2

### Generate random DAG
beta.min <- 0.5
beta.max <- 2
edge.pr <- 2 * ss / (pp - 1)
g <- pcalg::randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!

### Generate random data
pi <- sample(1:pp)
X <- pcalg::rmvDAG(n = nn, dag = g, errDist = "normal")
X <- X[, pi] ## permute the columns to randomize node ordering

test_that("Testing default behaviour of ccdr.run", {
    final <- ccdr.run(data = X, lambdas.length = 20)

    expect_is(final, "list")

    ### Check output types
    for(i in seq_along(final)){
        expect_is(final[[i]], "ccdrFit")
    }

    ### Check consistency of nedge
    for(i in seq_along(final)){
        matrix.nedge <- Matrix::nnzero(get.adjacency.matrix(final[[i]]$edges))
        edgeL.nedge <- num.edges(final[[i]]$edges)
        expect_equal(final[[i]]$nedge, edgeL.nedge, matrix.nedge)
    }
})

test_that("Testing ccdr.run with manual settings", {
    lambdas <- generate.lambdas(lambda.max = sqrt(nn), lambdas.ratio = 0.1, lambdas.length = 20, scale = "linear")
    final <- ccdr.run(data = X, lambdas = lambdas, alpha = 3, max.iters = 10, verbose = FALSE)

    expect_is(final, "list")

    ### Check output types
    for(i in seq_along(final)){
        expect_is(final[[i]], "ccdrFit")
    }

    ### Check consistency of nedge
    for(i in seq_along(final)){
        matrix.nedge <- Matrix::nnzero(get.adjacency.matrix(final[[i]]$edges))
        edgeL.nedge <- num.edges(final[[i]]$edges)
        expect_equal(final[[i]]$nedge, edgeL.nedge, matrix.nedge)
    }
})

### OLD TEST CODE
# source('~/Dropbox/PhD Research/Programming Projects/bncompare_dev/bncompare/R/bncompare-generate.R')
# # depends on
# source('~/Dropbox/PhD Research/Programming Projects/bncompare_dev/bncompare/R/s3-trueGraph.R')
# # depends on
# source('~/Dropbox/PhD Research/Programming Projects/bncompare_dev/bncompare/R/bncompare-utils.R')
# # depends on
# source('~/Dropbox/PhD Research/Programming Projects/bncompare_dev/bncompare/R/bncompare-functions.R')
#
# R_DEBUG_ON <<- FALSE # This is just a hack since we aren't using load_bncompare.R to load the package
#
# ### Generate some random data
# pp <- 10
# nn <- 100
# ss <- 2
#
# g <- generate_ordered_dag(pp, ss)
# d <- generate_data(g, nn)
# final <- ccdr.run(data = d$dat, lambdas.length = 20, alpha = 3, verbose = FALSE)
#
# test_that("Testing default behaviour of ccdr.gridR", {
#     expect_is(final, "list")
#
#     ### Check output types
#     for(i in seq_along(final)){
#         expect_is(final[[i]], "ccdrFit")
#     }
#
#     ### Check consistency of nedge
#     for(i in seq_along(final)){
#         matrix.nedge <- sum(as.matrix(final[[i]]$sbm) != 0)
#         sbm.nedge <- .num_edges(final[[i]]$sbm)
#         expect_equal(final[[i]]$nedge, sbm.nedge, matrix.nedge)
#     }
# })
