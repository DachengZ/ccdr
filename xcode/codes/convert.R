### functions to convert between matrix, edgeList, graph, and ccdrFit
### to graph
ccdrFit2graph <- function(cF) {
  cF.edL <- cF$edges
  cF.V <- as.character(1:cF$pp)
  names(cF.edL) <- cF.V
  return(graphNEL(nodes = cF.V, edgeL = cF.edL, edgemode = 'directed'))  
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
  ## o is the new order of nodes
  ## i is at position q[i] in vector o
  ## so edge i->j to q[i] -> q[j]
  q <- order(o)
  nodes0 <- nodes(g)
  edges0 <- edgeL(g)

  nn <- length(nodes0)
  ## nodes1 <- nodes0[o]
  
  edges1 <- vector("list", length = nn)
  names(edges1) <- nodes0
  for(i in 1:nn) {
    edges1[[q[i]]]$edges <- q[edges0[[i]]$edges]
    }
  return(graphNEL(nodes = nodes0, edgeL = edges1, edgemode = 'directed'))
}
