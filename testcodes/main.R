# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
# if(devtools::find_rtools()) devtools::install_github("DachengZ/ccdr")
# library(ccdri) ## manually build this package and load

setwd("~/ccdr/testcodes") ## change if necessary

# install.packages("bnlearn")
library(bnlearn)

# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
library(graph)

# install.packages("pcalg")
# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("RBGL")
# library(RBGL)
library(pcalg)

library(igraph)

source("rmvDAG.R") ## generate random data with optional intervention
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`
source("maintest.R") ## main test

### Set up the model parameters
nn <- 100                # How many samples to draw? ## this will be overridden if vfix specified
pp <- 50             # How many nodes in the DAG?
num.edges <- 50       # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 1
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
mm <- wgtMatrix(g, FALSE)

vfix <- rep(sample(1:pp), 5) # nodes to be fixed later ## this will override `nn`
N <- 1 # number of tests
test <- maintest(g, vfix, N = N)
# this will return average performance of repeated runs (now seems unnecessary)
colMeans(test$metric)


# nodes <- c(4, 40, 13, 39, 47) # change this to whatever nodes we want to add intervention
# vfix.rev <- rep(nodes[sample(length(nodes))], 50)
# test.rev <- maintest(g, vfix.rev, N = N, gamma = 2, originaldata = test$data, originalvfix = test$vfix)
