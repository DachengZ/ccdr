swaptest <- function(g, vfix = NULL, nn = NULL, N = 50, originaldata = NULL, originalvfix = NULL) {
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
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR","SHD")
    edges <- matrix(0, pp, pp) ## to count how many times each edges is estimated

    dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix)
    dat <- rbind(dat, originaldata) ## todo: check if sizes match
    vfix <- c(vfix, originalvfix)

    for(testi in 1:N) {
        ### Generate random data according to this DAG using the method rmvDAG

        o <- sample(1:pp)
        o1 <- c(o, pp + 1)
        q <- order(o) # o[q] == q[o] == 1:pp
        q1 <- order(o1)
        dat1 <- dat[, o]
        vfix1 <- q1[vfix] # permute the columns to randomize node ordering
        # permutations actually affects estimates

        ### Run the algorithm
        ccdr.path <- ccdr.run(data = dat1, intervention = vfix1, lambdas.length = 20, alpha = 10, verbose = FALSE)
        print(ccdr.path) # print some messages to check the speed of this algorithm

        ### Find the "best" one with the smallest SHD
        g1 <- permutenodes(g, o) ## currently this will lose weight data
        graph.path <- lapply(ccdr.path, ccdrFit2graph)
        compare.path <- sapply(graph.path, compare.graph, g1)
        shd.val <- compare.path[7, ]
        # lcp <- length(ccdr.path)
        # de <- rep(0, lcp)
        # for(i in 1:lcp) de[i] <- ccdr.path[[i]]$nedge
        # dr <- diff(shd.val) / diff(de)
        # z <- min(which(dr >= 0.4, arr.ind = TRUE))
        print(paste0(shd.val))
        z <- which(shd.val == min(shd.val))
        z <- z[length(z)]
        print(z)
        graph.shd <- permutenodes(graph.path[[z]], q)
        metric[testi, ] <- compare.path[, z]
        edges <- edges + wgtMatrix(graph.shd, transpose = FALSE)
    }
    ## includes the samples in the last test
    return(list(edges = edges, metric = metric, samples = nn, tests = N, data = dat, vfix = vfix))
}
