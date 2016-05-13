source("rmvDAG_fix.R") ## generate random data with points allowed to be fixed (in a specified range)
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`

maintest <- function(g, vfix = NULL, nn = NULL, N = 50, gamma = 2, originaldata = NULL, originalvfix = NULL) {
    ## g: input DAG
    ## vfix: nodes to be fixed in each run. length(vfix) == nn
    ## nn: samples to draw; only specify nn when vfix is NULL (no intervention)
    ## N: number of tests
    ## originaldata: previous observations/experiments data

    ## todo: distinguish nn and N from input?
    ## todo: make originaldata a list of data?

    pp <- numNodes(g)
    if(!is.null(vfix)) {
        nn <- length(vfix) ## todo: check error if both vfix and nn are NULL
    } else {
        if(is.null(nn)) stop("vfix and nn cannot both be NULL!") else vfix <- rep(pp + 1, nn)
    }
    if(!is.null(originaldata) & is.null(originalvfix)) originalvfix <- rep(pp + 1, nrow(originaldata))
    metric <- matrix(0, N, 7) ## to record performance metrics
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR", "SHD")
    edges <- matrix(0, pp, pp) ## to count how many times each edges is estimated
    time <- rep(0, N)

    for(testi in 1:N) {
        ### Generate random data according to this DAG using the method rmvDAG

        dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix)
        dat <- rbind(dat, originaldata) ## todo: check if sizes match

        o <- sample(1:pp)
        o1 <- c(o, pp + 1)
        q <- order(o) # o[q] == q[o] == 1:pp
        q1 <- order(o1)
        dat1 <- dat[, o] # permute the columns to randomize node ordering
        vfix1 <- q1[c(vfix, originalvfix)]
        # permutations actually affects estimates

        ### Run the algorithm
        ccdr.path <- ccdr.run(data = dat1, intervention = vfix1, gamma = gamma, lambdas.length = 20, alpha = 10, verbose = FALSE)
        print(ccdr.path) # print some messages to check the speed of this algorithm
        time[testi] <- sum(sapply(ccdr.path, getElement, "time"))

        ### Find the "best" one with the smallest SHD
        ### Or until SHD grows too fast
        g1 <- permutenodes(g, o)
        graph.path <- lapply(ccdr.path, ccdrFit2graph)
        compare.path <- sapply(graph.path, compare.graph, g1)
        shd.val <- compare.path[7, ]
        z <- which(shd.val == min(shd.val))
        z <- z[length(z)]
        lcp <- length(ccdr.path)
        if(z < lcp) {
            de <- rep(0, lcp - z + 1)
            for(i in 0:(lcp - z)) de[i + 1] <- ccdr.path[[z + i]]$nedge
            dr <- diff(shd.val[z:lcp]) / diff(de)
            z <- z + min(which(dr > 0.3, arr.ind = TRUE)) - 1
        }
        print(paste0(shd.val))
        print(z)
        graph.shd <- permutenodes(graph.path[[z]], q)
        metric[testi, ] <- compare.path[, z]
        edges <- edges + wgtMatrix(graph.shd, transpose = FALSE)
    }
    ## includes the samples in the last test
    return(list(edges = edges, metric = metric, time = sum(time), samples = nn, tests = N, data = dat, vfix = vfix))
}

summarynewtest <- function(test1, edges) {
    ## test1: new test (eg. with nodes fixed)
    ## edges: specify edges to be summarized; or based on vfix
    edges1 <- test1$edges
    return(cbind(from = edges[, 1], to = edges[, 2], est = edges1[edges[, 1:2, drop = FALSE]], rev = t(edges1)[edges[, 1:2, drop = FALSE]], old = edges[, 3], oldrev = edges[, 4]))
}

redge <- function(test) {
    edges <- test$edges
    redge0 <- which(t(edges) > 2 * edges & edges >= test$tests / 10, arr.ind = T)
    redge0 <- redge0[redge0[, 1] < redge0[, 2], ]
    redge0 <- matrix(redge0, ncol = 2)
    return(cbind(from = redge0[, 1], to = redge0[, 2], true = edges[redge0], rev = t(edges)[redge0]))
}

tedge <- function(test) {
    edges <- test$edges
    tedge0 <- which(edges > 2 * t(edges) & t(edges) >= test$tests / 10, arr.ind = T)
    tedge0 <- tedge0[tedge0[, 1] < tedge0[, 2], ]
    tedge0 <- matrix(tedge0, ncol = 2)
    return(cbind(from = tedge0[, 1], to = tedge0[, 2], true = edges[tedge0], rev = t(edges)[tedge0]))
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
