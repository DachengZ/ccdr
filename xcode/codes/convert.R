### functions to convert between matrix, edgeList, graph, and ccdrFit
### to graph
ccdrFit2graph <- function(cF) {
    ## edgelist of cF: children, then parent
    ## edgelist of graphNEL: parent, then children
    pp <- cF$pp
    cF.edL <- cF$edges
    edL <- vector("list", pp)
    for(i in 1:pp) edL[[i]] <- list(edges = NULL)
    for(i in 1:pp) {
        if(length(cF.edL[[i]]) > 0) {
            for(j in cF.edL[[i]]) edL[[j]]$edges <- c(edL[[j]]$edges, i)
        }
    }
    V <- as.character(1:pp)
    names(edL) <- V
    return(graphNEL(nodes = V, edgeL = edL, edgemode = 'directed'))
}

## to matrix
## m_{ij} = 1 <--> edge i->j exists
## or t(get.adjacency.matrix()) for a sparse matrix
## or wgtMatrix(, FALSE)
ccdrFit2matrix <- function(cF) {
    return(t(as(ccdrFit2graph(cF), "matrix") != 0))
}

## adjacency matrix to graph
## m_{ij} = 1 <--> edge i->j exists
adj2graph <- function(m) {
    return(as(graphAM(adj = t(m), edgemode = 'directed'), "graphNEL"))
}

## permute nodes (and edges) of a graph
## todo: add weight
permutenodes <- function(g, o) {
    ## so that when plot, they are the same
    ## but look at the matrix, the rows and columns are permuted
    ## o is the new order of nodes
    ## i is at position q[i] in vector o
    ## so edge i->j in the new graph is q[i] -> q[j] in the original graph
    q <- order(o)
    nodes0 <- nodes(g)
    edges0 <- edgeL(g)

    nn <- length(nodes0)
    nodes1 <- nodes0[o]

    edges1 <- vector("list", length = nn)
    names(edges1) <- nodes1
    for(i in 1:nn) {
        edges1[[q[i]]]$edges <- q[edges0[[i]]$edges]
    }
    return(graphNEL(nodes = nodes1, edgeL = edges1, edgemode = 'directed'))
}
