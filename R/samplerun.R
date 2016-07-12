# if(devtools::find_rtools()) devtools::install_github("DachengZ/ccdr")
# library(ccdri) ## manually build this package and load

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

### Set up the model parameters
nn <- 100               # How many samples to draw? ## this will be overridden if vfix specified
pp <- 50                # How many nodes in the DAG?
num.edges <- 50         # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 1
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
mm <- wgtMatrix(g, FALSE)

vfix <- rep(sample(1:pp), 5) # nodes to be fixed later ## this will override `nn`
## set vfix <- rep(pp + 1, nn) for no intervention
if(!is.null(vfix)) nn <- length(vfix)

dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix)

o <- sample(1:pp) ## to permute rows
o1 <- c(o, pp + 1)
q <- order(o) ## o[q] == q[o] == 1:pp
q1 <- order(o1)
dat1 <- dat[, o] # permute the columns to randomize node ordering
## permutations actually affects estimates

### Run the algorithm
ccdr.path <- ccdr.run(data = dat1, intervention = vfix, lambdas.length = 20, verbose = FALSE)
print(ccdr.path) # print solution path

### Find the "best" one with the smallest SHD
### Or until SHD grows too fast
g1 <- permutenodes(g, o)
graph.path <- lapply(ccdr.path, ccdrFit2graph)
compare.path <- sapply(graph.path, compare.graph, g1) ## to evaluate metrics such as TP, FP, TPR, SHD, etc.
shd.val <- compare.path[7, ]
z <- which(shd.val == min(shd.val))
z <- z[length(z)]
print(paste0(shd.val))

graph.shd <- permutenodes(graph.path[[z]], q) ## this is the "best" estimate
