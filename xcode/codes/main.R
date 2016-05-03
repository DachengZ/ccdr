# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
library(ccdr)

setwd("~/ccdr/xcode/codes")

# install.packages("bnlearn")
library(bnlearn)

# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
library(graph)

# install.packages("pcalg)
# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("RBGL")
# library(RBGL)
library(pcalg)

#library(igraph)

source("rmvDAG_fix.R") ## generate random data with points allowed to be fixed (in a specified range)
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`
source("maintest.R") ## main test
source("swaptest.R")

### Set up the model parameters
nn <- 100                # How many samples to draw?
pp <- 50              # How many nodes in the DAG?
num.edges <- 100       # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 2
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
mm <- wgtMatrix(g, FALSE)

vfix <- c() # nodes to be fixed later
N <- 50 # number of tests
test <- swaptest(g, vfix, nn, N)

## intervention on all points
vfix <- rep(sample(pp), 5)
test1 <- maintest(g, vfix, N = N)


## focus on reversed edges
## get edges with fewer true estimates and far more reversed estimates
redge0 <- redge(test)
redge0 ## why nothing
tedge0 <- tedge(test)
tedge0 ## why nothing

# redge1 <- redge(test1)
# redge1 ## why things

summarynewtest(test1, redge0)
# summarynewtest(test, redge1)

## fix -outgoing- nodes from reversed edges
revnodes <- c(redge0[, 1], redge0[, 2]) ## unique?

revnodes <- c(17, 26, 30)
for(j in revnodes) {
    revnodes <- c(revnodes, as.integer(inEdges(as.character(j), g)[[1]]))
}
#revnodes <- c(10, 19, 34, 41, 14)
vfix.rev <- rep(revnodes[sample(length(revnodes))], 2 * nn)

# test.rev <- swaptest(g, vfix.rev, N = N)
test.rev <- swaptest(g, vfix.rev, N = N, originaldata = test$data)

## see if any changes
summarynewtest(test.rev, redge0)
summarynewtest(test.rev, tedge0)


## usually the edges are still in estimates, and some (most) estimates are still reversed
## this heavily depends on how we fix the nodes
## eg. what distribution? (runif)
## what range? (9 - 11)
## if the incoming nodes has more than one parents, do we fix other parent nodes?

### or fixing incoming nodes
### usually the edges no longer appear in any estimates
### so it is less likely that a causal relationship exists in this direction

## non-existent edge
## Case I: a->b->c and we find a->c(c->a)
## Case II: a->b, c->b and we find a->c(c->a)
nedge0 <- nedge(test)
vfix.new.i <- unique(nedge0[, 1])
test.new.i <- maintest(g, vfix.new.i, nn, N)
summarynewtest(test, test.new.i)
