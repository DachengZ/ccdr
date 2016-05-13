library(gdata)
# praf, pmek, plcg, PIP2, PIP3, p44/42, pakts473, PKA, PKC, P38, pjnk
##
m <- matrix(0, 11, 11)
index <- cbind(c(1, 2, 3, 3, 4, 5, 5, 5, 6, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9), c(2, 6, 4, 9, 9, 3, 4, 7, 7, 1, 2, 6, 7, 10, 11, 1, 2, 8, 10, 11))
m[index] <- 1
g <- as(graphAM(adj = m, edgemode = 'directed'), "graphNEL")
nodename <- c("praf", "pmek", "plcg", "PIP2", "PIP3", "p44/42", "pakts473", "PKA", "PKC", "P38", "pjnk")
nodes(g) <- nodename
vfix <- c()

# data.null <- read.xls("../data/1.\ cd3cd28.xls")
# vfix <- c(vfix, rep(12, dim(data.null)[1]))

# data.null1 <- read.xls("../data/2.\ cd3cd28icam2.xls")
# vfix <- c(vfix, rep(12, dim(data.null1)[1])

data.pka <- read.xls("../data/9.\ b2camp.xls")
vfix <- c(vfix, rep(8, dim(data.pka)[1]))

data.akt <- read.xls("../data/3.\ cd3cd28+aktinhib.xls")
vfix <- c(vfix, rep(7, dim(data.akt)[1]))

data.akt1 <- read.xls("../data/10.\ cd3cd28icam2+aktinhib.xls")
vfix <- c(vfix, rep(7, dim(data.akt1)[1]))

data.mek <- read.xls("../data/6.\ cd3cd28+u0126.xls")
vfix <- c(vfix, rep(6, dim(data.mek)[1])) ## 2 or 6?

data.mek1 <- read.xls("../data/13.\ cd3cd28icam2+u0126.xls")
vfix <- c(vfix, rep(6, dim(data.mek1)[1])) ## 2 or 6

data.pkc <- read.xls("../data/8.\ pma.xls")
colnames(data.pkc) <- colnames(data.pka)
vfix <- c(vfix, rep(9, dim(data.pkc)[1]))

data.pkc1 <- read.xls("../data/4.\ cd3cd28+g0076.xls")
vfix <- c(vfix, rep(9, dim(data.pkc1)[1]))

data.pkc2 <- read.xls("../data/11.\ cd3cd28icam2+g0076.xls")
vfix <- c(vfix, rep(9, dim(data.pkc2)[1]))

data.pip2 <- read.xls("../data/5.\ cd3cd28+psitect.xls")
vfix <- c(vfix, rep(4, dim(data.pip2)[1]))

data.pip21 <- read.xls("../data/12.\ cd3cd28icam2+psit.xls")
vfix <- c(vfix, rep(4, dim(data.pip21)[1]))

data.akt2 <- read.xls("../data/7.\ cd3cd28+ly.xls")
vfix <- c(vfix, rep(7, dim(data.akt2)[1]))

data.akt3 <- read.xls("../data/14.\ cd3cd28icam2+ly.xls")
vfix <- c(vfix, rep(7, dim(data.akt3)[1]))

data <- rbind(#data.null, data.null1,
              data.pka, data.akt, data.akt1,
              data.mek, data.mek1,
              data.pkc,
              data.pkc1, data.pkc2,
              data.pip2, data.pip21,
              data.akt2, data.akt3
              )
data <- data.matrix(data)

data <- log(data)
vfix <- as.integer(vfix)

pp <- 11
N <- 10
metric <- matrix(0, N, 7) ## to record performance metrics
colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR")
edges <- matrix(0, pp, pp) ## to count how many times each edges is estimated
testi <- 1
for(testi in 1:N) {
    o <- sample(1:pp)
    o1 <- c(o, pp + 1)
    q <- order(o) # o[q] == q[o] == 1:pp
    q1 <- order(o1)
    dat1 <- data[, o] # permute the columns to randomize node ordering
    vfix1 <- q1[vfix]

    ccdr.path <- ccdr.run(data = dat1, intervention = vfix1, lambdas.length = 20, gamma = -1, alpha = 10, verbose = FALSE)
    print(ccdr.path)

    g1 <- permutenodes(g, o)
    graph.path <- lapply(ccdr.path, ccdrFit2graph)
    compare.path <- sapply(graph.path, compare.graph, g1)
    shd.val <- compare.path[7, ]
    print(paste0(shd.val))

    npath <- length(ccdr.path)
    nedge <- rep(0, npath)
    for(i in 1:npath) nedge[i] <- ccdr.path[[i]]$nedge
    z0 <- which(abs(nedge - 20) == min(abs(nedge - 20)), arr.ind = T)
    z1 <- which(shd.val[z0] == min(shd.val[z0]), arr.ind = T)
    z <- z0[length(z1)]
    print(shd.val[z])

    graph.shd <- permutenodes(graph.path[[z]], q)
    nodes(graph.shd) <- nodename
    metric[testi, ] <- compare.graph(graph.shd, g)
    edges <- edges + wgtMatrix(graph.shd, transpose = FALSE)

}


## to plot
ig <- igraph.from.graphNEL(g)
ig.c <- layout_in_circle(ig)
par(mfrow = c(1, 2))
plot(ig, layout = ig.c, vertex.color = "white", edge.color = "black", edge.arrow.size = 0.7)
plot(igraph.from.graphNEL(graph.shd), layout = ig.c, vertex.color = "white", edge.color = "black", edge.arrow.size = 0.7)
