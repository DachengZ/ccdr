swaptest <- function(g, vfix = NULL, nn = NULL, N = 50, originaldata = NULL) {
    ## g: input DAG
    ## vfix: nodes to be fixed in each run. length(vfix) == nn
    ## nn: samples to draw; only specify nn when vfix is NULL (no intervention)
    ## N: number of tests
    ## originaldata: previous observations/experiments data

    ## todo: distinguish nn and N from input?
    ## todo: make originaldata a list of data?

    pp <- length(g@nodes)
    if(!is.null(vfix)) nn <- length(vfix) ## todo: check error if both vfix and nn are NULL
    else vfix <- rep(pp + 1, nn)
    metric <- matrix(0, N, 6) ## to record performance metrics
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR")
    edges <- matrix(0, pp, pp) ## to count how many times each edges is estimated

    dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix)
    dat <- rbind(dat, originaldata) ## todo: check if sizes match
    vfix <- dat[, pp + 1]

    for(testi in 1:N) {
        ### Generate random data according to this DAG using the method rmvDAG

        o <- sample(1:pp)
        o1 <- c(o, pp + 1)
        q <- order(o) # o[q] == q[o] == 1:pp
        q1 <- order(o1)
        dat1 <- cbind(dat[, o], q1[vfix]) # permute the columns to randomize node ordering
        # permutations actually affects estimates

        ### Run the algorithm
        ccdr.path <- ccdr.run(data = dat1, lambdas.length = 20, alpha = 10, verbose = FALSE)
        print(ccdr.path) # print some messages to check the speed of this algorithm

        ### Find the "best" one with the smallest SHD
        g1 <- permutenodes(g, o) ## currently this will lose weight data
        graph.path <- lapply(ccdr.path, ccdrFit2graph)
        shd.val <- sapply(graph.path, pcalg::shd, g1)
        lcp <- length(ccdr.path)
        de <- rep(0, lcp)
        for(i in 1:lcp) de[i] <- ccdr.path[[i]]$nedge
        dr <- diff(shd.val) / diff(de)
        z <- min(which(dr > 0.5, arr.ind = TRUE))
        print(paste0(shd.val))
        print(z)
        # z <- which(shd.val == min(shd.val))
        # z <- z[length(z)]
        graph.shd <- permutenodes(graph.path[[z]], q)
        metric[testi, ] <- compare.graph(graph.shd, g)
        edges <- edges + wgtMatrix(graph.shd, transpose = FALSE)
    }
    ## includes the samples in the last test
    return(list(edges = edges, metric = metric, samples = nn, tests = N, data = dat))
}
