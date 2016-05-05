# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
# library(ccdr) ## manually build this package and load

setwd("~/ccdr/xcode/codes") ## change if necessary

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
## source("swaptest.R") ## outdated method if in each test we just swap columns on the same sample

### Set up the model parameters
nn <- 50                # How many samples to draw?
pp <- 100              # How many nodes in the DAG?
num.edges <- 50       # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 2
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
mm <- wgtMatrix(g, FALSE)

vfix <- c() # nodes to be fixed later
N <- 50 # number of tests
test <- maintest(g, vfix, nn, N)
colMeans(test$metric)
apply(test$metric, 2, sd)

## intervention on all points
vfix <- sample(1:(pp+1), nn, replace = T)
# test1 <- maintest(g, vfix, N = N)
test1 <- maintest(g, vfix, N = N, originaldata = test$data)

## focus on reversed edges
## get edges with fewer true estimates and many more reversed estimates
redge0 <- redge(test)
redge0
## and true edges
tedge0 <- tedge(test)
tedge0

summarynewtest(test1, redge0)
summarynewtest(test1, tedge0)
# summarynewtest(test, redge1)

## fix -outgoing- nodes from reversed edges
## revnodes <- c(redge0[, 1], redge0[, 2]) ## unique?
## revnodes <- c(15, 19)
## for(j in revnodes) {
##     revnodes <- c(revnodes, as.integer(inEdges(as.character(j), g)[[1]]))
## }
revnodes <- c(15, 30, 2, 49)
vfix.rev <- rep(revnodes[sample(length(revnodes))], nn)
test.rev <- maintest(g, vfix.rev, N = N)
# test.rev <- maintest(g, vfix.rev, N = N, originaldata = test$data)

colMeans(test.rev$metric)
apply(test.rev$metric, 2, sd)
## see if any changes
summarynewtest(test.rev, redge0)
summarynewtest(test.rev, tedge0)
