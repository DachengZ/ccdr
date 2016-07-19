# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
# if(devtools::find_rtools()) devtools::install_github("DachengZ/ccdr")
library(ccdr)
library(ccdri) ## manually build this package and load

setwd("~/ccdr/testcodes") ## change if necessary

library(bnlearn)
library(graph)
library(RBGL)
library(pcalg)
library(igraph)

source("rmvDAG.R") ## generate random data with optional intervention
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`
source("maintest.R") ## main test

### Test function wrapper
ccdrtest <- function(nn, pp, num.edges, itvtimes = NA, vfix = NULL, gamma = 2, lambdas.length = 20) {
    # Goal: Based on a same graph, we generata two sets of data,
    #       one observational, one with intervention
    #       And test CCDr, CCDri, PC on them.

    ## Parameters for the DAG
    # nn: number of samples to draw for observational data
    # pp: number of many nodes in the DAG
    # num.edges: number of *expected* edges in the DAG
    ss <- num.edges / pp # the expected number of parents *per node*

    ## To record performance
    metric <- matrix(0, 6, 8) ## to record performance metrics

    ### Generate a random DAG using the pcalg method randomDAG
    beta.min <- 0.5
    beta.max <- 1
    edge.pr <- 2 * ss / (pp - 1)
    g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!

    ### Observational data
    vfix0 <- as.integer(rep(pp + 1, nn))
    dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix0)

    ## Permute columns
    o <- sample(1:pp)
    o1 <- c(o, pp + 1)
    q <- order(o) ## o[q] == q[o] == 1:pp
    q1 <- order(o1)
    dat1 <- dat[, o] # permute the columns to randomize node ordering
    g1 <- permutenodes(g, o)
    ## permutations actually affects estimates

    ## CCDr
    ptm <- proc.time()
    ccdr.path <- ccdr::ccdr.run(data = dat1, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdr.path)
    graph.path <- lapply(ccdr.path, ccdrFit2graph)
    compare.path <- sapply(graph.path, compare.graph, g1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    # graph.shd <- permutenodes(graph.path[[z]], q)
    metric[1, 1:7] <- compare.path[, z]
    metric[1, 8] <- (proc.time() - ptm)[3]
    # metric[1, 8:10] <- (proc.time() - ptm)[1:3]

    ## CCDri
    ptm <- proc.time()
    ccdri.path <- ccdri::ccdr.run(data = dat1, intervention = vfix0, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdri.path)
    graph.path <- lapply(ccdri.path, ccdrFit2graph)
    compare.path <- sapply(graph.path, compare.graph, g1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    # graph.shd <- permutenodes(graph.path[[z]], q)
    metric[2, 1:7] <- compare.path[, z]
    metric[2, 8] <- (proc.time() - ptm)[3]
    # metric[2, 8:10] <- (proc.time() - ptm)[1:3]

    ## MMHC from `bnlearn`
    # PC from `pcalg` gives bi-directed DAG # Why?
    ptm <- proc.time()
    # pc.fit <- pc(suffStat = list(C = cor(dat1), n = nn1), indepTest = gaussCItest, alpha = 0.01, p = pp, verbose = F)
    mmhc.fit <- mmhc(as.data.frame(dat1))
    mmhc.graph <- as.graphNEL(mmhc.fit)
    metric[3, 1:7] <- compare.graph(mmhc.graph, g1)
    metric[3, 8] <- (proc.time() - ptm)[3]
    # metric[3, 8:10] <- (proc.time() - ptm)[1:3]

    true.edges <- metric[1, 2] - metric[1, 4] + metric[1, 7]

    ### Intervention data
    if(is.null(vfix)) {
        if(is.na(itvtimes)) stop("Please provide either vfix or itvtimes")
        else vfix <- rep(sample(pp), itvtimes)
    } else itvtimes <- paste("~", length(vfix) / pp, sep = "")
    # it seems that, here itvtimes can be non-integer;
    # R will replace itvtimes by floor(itvtimes)
    nn1 <- length(vfix)
    dat <- rmvDAG.fix(n = nn1, dag = g, vfix = vfix)
    dat1 <- dat[, o] # permute the columns to randomize node ordering
    vfix1 <- q1[vfix]
    g1 <- permutenodes(g, o)
    ## permutations actually affects estimates

    ## CCDr
    ptm <- proc.time()
    ccdr.path <- ccdr::ccdr.run(data = dat1, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdr.path)
    graph.path <- lapply(ccdr.path, ccdrFit2graph)
    compare.path <- sapply(graph.path, compare.graph, g1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    # graph.shd <- permutenodes(graph.path[[z]], q)
    metric[4, 1:7] <- compare.path[, z]
    metric[4, 8] <- (proc.time() - ptm)[3]
    # metric[4, 8:10] <- (proc.time() - ptm)[1:3]

    ## CCDri
    ptm <- proc.time()
    ccdri.path <- ccdri::ccdr.run(data = dat1, intervention = vfix1, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdri.path)
    graph.path <- lapply(ccdri.path, ccdrFit2graph)
    compare.path <- sapply(graph.path, compare.graph, g1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    # graph.shd <- permutenodes(graph.path[[z]], q)
    metric[5, 1:7] <- compare.path[, z]
    metric[5, 8] <- (proc.time() - ptm)[3]
    # metric[5, 8:10] <- (proc.time() - ptm)[1:3]

    ## MMHC from `bnlearn`
    # PC from `pcalg` gives bi-directed DAG # Why?
    ptm <- proc.time()
    # pc.fit <- pc(suffStat = list(C = cor(dat1), n = nn1), indepTest = gaussCItest, alpha = 0.01, p = pp, verbose = F)
    mmhc.fit <- mmhc(as.data.frame(dat1))
    mmhc.graph <- as.graphNEL(mmhc.fit)
    metric[6, 1:7] <- compare.graph(mmhc.graph, g1)
    metric[6, 8] <- (proc.time() - ptm)[3]
    # metric[6, 8:10] <- (proc.time() - ptm)[1:3]

    namelist <- rep(c("CCDr", "CCDri", "MMHC"), 2)
    nnlist <- c(rep(nn, 3), rep(length(vfix), 3))
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR", "SHD", "elapsed")
    itvlist <- c(rep("obs", 3), rep("itv", 3))
    itvtimes <- c(rep(NA, 3), rep(itvtimes, 3)) # use 0 instead of NA?
    # pplist <- as.integer(rep(pp, 6))
    # nelist <- as.integer(rep(num.edges, 6))
    return(data.frame(name = namelist, nn = nnlist, pp, num.edges,
                      itv = itvlist, itvtimes = itvtimes, true.edges,
                      metric))
}

m <- vector("list", 10)

m[[1]] <- ccdrtest(nn = 20, pp = 10, num.edges = 10, itvtimes = 2)
m[[2]] <- ccdrtest(nn = 50, pp = 10, num.edges = 10, itvtimes = 5)
m[[3]] <- ccdrtest(nn = 50, pp = 20, num.edges = 20, itvtimes = 2)
m[[4]] <- ccdrtest(nn = 100, pp = 20, num.edges = 20, itvtimes = 5)
m[[5]] <- ccdrtest(nn = 100, pp = 50, num.edges = 20, itvtimes = 2)
m[[6]] <- ccdrtest(nn = 100, pp = 50, num.edges = 50, itvtimes = 2)
m[[7]] <- ccdrtest(nn = 200, pp = 100, num.edges = 50, itvtimes = 2)
m[[8]] <- ccdrtest(nn = 200, pp = 100, num.edges = 100, itvtimes = 2)
m[[9]] <- ccdrtest(nn = 500, pp = 200, num.edges = 100, itvtimes = 2)
m[[10]] <- ccdrtest(nn = 500, pp = 200, num.edges = 200, itvtimes = 2)

# CCDri performs well when intervention times is the same for each node (?)


# The following setting takes over 5 minutes on intervention data.
# m[[10]] <- ccdrtest(nn = 500, pp = 500, num.edges = 500, vfix = sample(500, 1000, T))

getm <- function(row, column) {
    return(sapply(m, "[", row, column))
}

## Check one simulation:
# m[[5]]

pdf(file = "obs.pdf")
par(oma = rep(0, 4), mar = c(4, 4, 2, 1))
plot(getm(1, 14), type = "p", pch = 1, main = "Timing Observational Data", xlab = "test index", ylab = "Time elapsed (seconds)", ylim = c(0, 5))
points(getm(1, 14), type = "l", lty = 1)
points(getm(2, 14), type = "p", pch = 20)
points(getm(2, 14), type = "l", lty = 2)
points(getm(3, 14), type = "p", pch = 24)
points(getm(3, 14), type = "l", lty = 3)
legend(1, 5, legend = c("CCDr", "CCDri", "MMHC"), lty = c(1, 2, 3), pch = c(1, 20, 24))
dev.off()

pdf(file = "itv.pdf")
par(oma = rep(0, 4), mar = c(4, 4, 2, 1))
plot(getm(4, 14), type = "p", pch = 1, main = "Timing Interventional Data", xlab = "test index", ylab = "Time elapsed (seconds)", ylim = c(0, 5))
points(getm(4, 14), type = "l", lty = 1)
points(getm(5, 14), type = "p", pch = 20)
points(getm(5, 14), type = "l", lty = 2)
points(getm(6, 14), type = "p", pch = 24)
points(getm(6, 14), type = "l", lty = 3)
legend(1, 5, legend = c("CCDr", "CCDri", "MMHC"), lty = c(1, 2, 3), pch = c(1, 20, 24))
dev.off()
