# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
library(ccdr)

setwd("~/Dropbox/codes")

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

### Set up the model parameters
nn <- 50                # How many samples to draw?
pp <- 50              # How many nodes in the DAG?
num.edges <- 25       # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 2
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
mm <- wgtMatrix(g, FALSE)

vfix <- c() # nodes to be fixed later
N <- 50 # number of tests
test <- maintest(g, vfix, nn)

## we found error
## intervention on all points
vfix <- rep(sample(pp), 5)
test1 <- maintest(g, vfix, 5 * pp)


## focus on reversed edges
## get edges with fewer true estimates and far more reversed estimates
redge0 <- redge(test)
redge0 ## why nothing
redge1 <- redge(test1)
redge1 ## why things

# summarynewtest(test1, redge0)
summarynewtest(test, redge1)

## fix -outgoing- nodes from reversed edges
# rev.o <- unique(redge0[, 1]) # cautious when some nodes are both incoming and outgoing
revnodes <- unique(c(redge1[, 1], redge1[, 2]))
vfix.rev <- rep(revnodes[sample(length(revnodes))], 10)

test.rev <-maintest(g, vfix.rev)
test.rev <- maintest(g, vfix.rev, originaldata = test$data)

## see if any changes
summarynewtest(test.rev, redge0)

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
