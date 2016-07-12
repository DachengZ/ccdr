### how to merge small graphs
makegraph <- function(p, e, K, beta.min = 0.5, beta.max = 1) {
    # p: number of total nodes
    # e: number of expected total edges
    # K: number of (expected) components
    # TODO: allow K to be a vector of the number of nodes in each component
    label <- sample(K, p, replace = T)
    pp <- rep(0, K)
    for(i in 1:K) pp[i] <- sum(label == i) ## or use table(a)
    ee <- e * pp / p ## expected number of edges in each component, can be decimal
    ss <- e / p # This is the expected number of parents *per node*

    ### Generate a random DAG using the pcalg method randomDAG
    nodesum <- 0
    glist <- vector(mode = "list", length = K)
    for(i in 1:K) {
        if(pp[i] > 1) {
            edge.pr <- 2 * ss / (pp[i] - 1)
            gi <- randomDAG(n = pp[i], prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!
            nodes(gi) <- as.character(nodesum + as.numeric(nodes(gi)))
            glist[[i]] <- graph_from_graphnel(gi)
        } else {
            gi <- graphNEL(nodes = as.character(nodesum + 1), edgemode = 'directed')
            glist[[i]] <- graph_from_graphnel(gi)
        }
        nodesum <- nodesum + pp[i]
    }
    gg <- do.call(union, glist)
    return(as_graphnel(gg))
}
