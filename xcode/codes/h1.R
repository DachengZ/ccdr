# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
# library(ccdr) ## manually build this package and load

setwd("~/ccdr/xcode/codes") ## change if necessary

## test on established graphs from `bnlearn` package
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

source("rmvDAG_fix.R") ## generate random data with points allowed to be fixed (in a specified range)
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`
source("maintest.R") ## main test
source("swaptest.R") ## outdated method if in each test we just swap columns on the same sample

# hailfinder, pathfinder, munin4/munin1, andes, win95pts
# andes underestimate but low FDR
# win95pts intervention has little effect
dataname <- "win95pts"
dataurl <- paste0("http://www.bnlearn.com/bnrepository/", dataname, "/", dataname, ".rda")
con <- url(dataurl)
load(con)
close(con)

### Convert the bnlearn data to a graph object
madj <- amat(bn)
g <- as(madj, "graphNEL")

pp <- numNodes(g)
nodename <- nodes(g)
nodes(g) <- as.character(1:pp)
### if not topologically sorted
### assume node are named by numbers 1:pp
### g1 <- permutenodes(g, as.integer(tsort(g)))

# change weight?
# default is 1. change to 0.5~1?
g <- changeweight(g, 0.5, 1)
vfix <- rep(sample(1:pp), 5)
# vfix <- sample(1:pp, 5, replace = T)
test <- maintest(g, vfix = vfix, N = 50)

## focus on reversed edges
## get edges with fewer true estimates and many more reversed estimates
redge0 <- redge(test)
redge0
## and true edges
tedge0 <- tedge(test)
tedge0

revnodes <- NULL
revnodes <- c(2, 3, 12, 5, 44, 48, 31, 46, 47, 49, 42, 54, 8, 56, 58, 39, 72) # change this to whatever nodes we want to add intervention
## if intervention only on the both ends of an edge has little effect,
## probably because there are other edges towards the receiving node
## try intervention on some other nodes that have edges towards the receving node as well
## or because the original data is "bad"
## how do we determine if intervention is "effective"?
## 50/0? 40/5? 30/10? 20/5?
vfix.rev <- rep(revnodes[sample(length(revnodes))], 20)
# test.rev0 <- maintest(g, vfix.rev, N = 50)
test.rev <- maintest(g, vfix.rev, N = 50, originaldata = test$data, originalvfix = test$vfix)

colMeans(test.rev$metric)
apply(test.rev$metric, 2, sd)
## see if any changes
summarynewtest(test.rev, redge0)
summarynewtest(test.rev, tedge0)

test.null <- maintest(g, n = 200, N = 50)
