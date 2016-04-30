source("rmvDAG_fix.R") ## generate random data with points allowed to be fixed (in a specified range)
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`

maintest <- function(g, vfix = NULL, nn = NULL, N = 50, originaldata = NULL) {
    ## g: input DAG
    ## vfix: nodes to be fixed in each run. length(vfix) == nn
    ## nn: samples to draw; only specify nn when vfix is NULL (no intervention)
    ## N: number of tests
    ## originaldata: previous observations/experiments data

    ## todo: distinguish nn and N from input?
    ## todo: make originaldata a list of data?

    if(!is.null(vfix)) nn <- length(vfix) ## todo: check error if both vfix and nn are NULL
                                          ## then only do test with the data we have in originaldata
    pp <- length(g@nodes)
    metric <- matrix(0, N, 6) ## to record performance metrics
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR")
    edges <- matrix(0, pp, pp) ## to count how many times each edges is estimated

    for(testi in 1:N) {
        ### Generate random data according to this DAG using the method rmvDAG

        dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix)
        dat <- rbind(dat, originaldata) ## todo: check if sizes match
        o <- sample(1:pp)
        q <- order(o) # o[q] == q[o] == 1:pp
        dat1 <- dat[, c(o, pp + 1)] # permute the columns to randomize node ordering
        # permutations actually affects estimates

        ### Run the algorithm
        ccdr.path <- ccdr.run(data = dat1, lambdas.length = 20, alpha = 10, verbose = TRUE)
        print(ccdr.path) # print some messages to check the speed of this algorithm

        ### Find the "best" one with the smallest SHD
        g1 <- permutenodes(g, o)
        graph.path <- lapply(ccdr.path, ccdrFit2graph)
        shd.val <- sapply(graph.path, pcalg::shd, g1)
        z <- which(shd.val == min(shd.val))
        z <- z[length(z)]
        graph.shd <- permutenodes(graph.path[[z]], q)

        #metric[testi, ] <- compare.graph(graph.shd, g)
        compare.graph(graph.shd, g)
        edges <- edges + wgtMatrix(graph.shd, transpose = FALSE)
    }
    return(list(edges = edges, metric = metric, samples = nn, tests = N, data = dat))
}

summarynewtest <- function(test1, edges) {
    ## test1: new test (eg. with nodes fixed)
    ## edges: specify edges to be summarized; or based on vfix
    edges1 <- test1$edges
    return(cbind(from = edges[, 1], to = edges[, 2], est = edges1[edges[, 1:2]], rev = t(edges1)[edges[, 1:2]], old = edges[, 3], oldrev = edges[, 4]))
}

summarynewtest.old <- function(test1, edges = NULL, test0 = NULL) {
    ## test0: original test
    ## test1: new test (eg. with nodes fixed)
    ## edges: specify edges to be summarized; or based on vfix
    edges1 <- test1$edges
    if(is.null(edges)) {
        vfix <- test1$vfix
        edges0 <- test0$edges
        ## first, see vfix as outgoing nodes, find their child nodes
        out <- which(edges0[vfix, , drop = FALSE] != 0, arr.ind = T)
        out[, 1] <- vfix[out[, 1]]
        out <- matrix(out[out[, 1] < out[, 2], ], ncol = 2)
        out1 <- cbind(from = out[, 1], to = out[, 2], est = edges1[out], rev = t(edges1)[out])

        ## then, see vfix as incoming nodes
        inc <- which(edges0[, vfix, drop = FALSE] != 0, arr.ind = T)
        inc[, 2] <- vfix[inc[, 2]]
        inc <- matrix(inc[inc[, 1] < inc[, 2], ], ncol = 2)
        inc1 <- cbind(from = inc[, 1], to = inc[, 2], est = edges1[inc], rev = t(edges1)[inc])

        return(rbind(out1, inc1))
    } else {
        return(cbind(from = edges[, 1], to = edges[, 2], est = edges1[edges[, 1:2]], rev = t(edges1)[edges[, 1:2]]))
    }
}

redge <- function(test) {
    edges <- test$edges
    redge0 <- which(t(edges) > 2 * edges & edges >= test$tests / 10, arr.ind = T)
    redge0 <- redge0[redge0[, 1] < redge0[, 2], ]
    redge0 <- matrix(redge0, ncol = 2)
    return(cbind(from = redge0[, 1], to = redge0[, 2], true = edges[redge0], rev = t(edges)[redge0]))
}

nedge <- function(test, b = NULL) {
    edges <- test$edges
    if(is.null(b)) b <- 0.4 * test$tests
    nedge0 <- which((mm == 0) & (edges + t(edges) != 0), arr.ind = T)
    ## in reality how do we get mm?
    nedge0 <- nedge0[nedge0[, 1] < nedge0[, 2] & edges[nedge0] + t(edges)[nedge0] >= b]
    nedge0 <- matrix(nedge0, ncol = 2)
    nedge0 <- cbind(from = nedge0[, 1], to = nedge0[, 2], est = edges[nedge0], rev = t(edges)[nedge0])
    return(nedge0)
}
