### scalefree

ig <- sample_pa(200)
g0 <- as_graphnel(ig)
g <- permutenodes(g0, as.integer(tsort(g0)))

pp <- numNodes(g)
vfix <- rep(sample(1:pp), 5)
test <- maintest(g, vfix = vfix, N = 10)

## focus on reversed edges
## get edges with fewer true estimates and many more reversed estimates
redge0 <- redge(test)
redge0

revnodes <- c(28, 32, 63, 71, 162, 179) # change this to whatever nodes we want to add intervention
## if intervention only on the both ends of an edge has little effect,
## probably because there are other edges towards the receiving node
## try intervention on some other nodes that have edges towards the receving node as well
## or because the original data is "bad"
## how do we determine if intervention is "effective"?
## 50/0? 40/5? 30/10? 20/5?
vfix.rev <- rep(revnodes[sample(length(revnodes))], 50)
# test.rev0 <- maintest(g, vfix.rev, N = 50)
test.rev <- maintest(g, vfix.rev, N = 20, originaldata = test$data, originalvfix = test$vfix)

colMeans(test.rev$metric)
apply(test.rev$metric, 2, sd)
## see if any changes
summarynewtest(test.rev, redge0)
