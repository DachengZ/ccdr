//
//  main.cpp
//  xcode
//
//  Created by Dacheng Zhang on 16/4/29.
//  Copyright © 2016年 Dacheng Zhang. All rights reserved.
//

//#include <iostream>
//
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    std::cout << "Hello, World!\n";
//    return 0;
//}

#include  <RInside.h>                   // for the embedded R via RInside
#include "../src/algorithm.h"

Rcpp::NumericMatrix createMatrix(const int n) {
    Rcpp::NumericMatrix M(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            M(i,j) = i*10 + j;
        }
    }
    return(M);
}

int main(int argc, char *argv[]) {
    
    RInside R(argc, argv);                      // create an embedded R instance
    
/*    const int mdim = 4;                         // let the matrices be 4 by 4; create, fill
    R["M"] = createMatrix(mdim);                // then assign data Matrix to R's 'M' var
    
    std::string str =
    "cat('Running ls()\n'); print(ls()); "
    "cat('Showing M\n'); print(M); "
    "cat('Showing colSums()\n'); Z <- colSums(M); print(Z); "
    "Z";                     // returns Z
    
    
    Rcpp::NumericVector v = R.parseEval(str);   // eval string, Z then assigned to num. vec
    
    for (int i=0; i< v.size(); i++) {           // show the result
        std::cout << "In C++ element " << i << " is " << v[i] << std::endl;
    }
*/
    std::string cmd =
    //"setwd(\"~/ccdr\")"
    "library(ccdr);"
    "library(bnlearn);"
    "library(graph);"
    "library(pcalg);"
    
    "setwd(\"~/ccdr/xcode/codes\");"
    
    "print(\"loading rmvDAG_fix\");"
    "source(\"rmvDAG_fix.R\");"

    "print(\"loading convert\");"
    "source(\"convert.R\");"

    "print(\"loading compare\");"
    "source(\"compare.R\");"

    "print(\"loading maintest\");"
    "source(\"maintest.R\");"

    //"print(\"loading ccdr-main-R\");"
    //"source(\"~/ccdr/R/ccdr-main-R.R\");"
    
    "nn <- 50;"
    "pp <- 50;"
    "num.edges <- 25;"
    "ss <- num.edges / pp;"
    "beta.min <- 0.5;"
    "beta.max <- 2;"
    "edge.pr <- 2 * ss / (pp - 1);"
    "g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max);"
    "mm <- wgtMatrix(g, FALSE);"
    "vfix <- c();"
    "N <- 50;"
    "pp <- length(g@nodes);"
    "dat <- rmvDAG.fix(n = nn, dag = g, vfix = vfix);"
    "o <- sample(1:pp);"
    "q <- order(o);"
    "dat1 <- dat[, c(o, pp + 1)];"
    "ccdr.path <- ccdr.run(data = dat1, lambdas.length = 20, alpha = 10, verbose = TRUE);"
    "print(ccdr.path);"
    "g1 <- permutenodes(g, o);"
    "graph.path <- lapply(ccdr.path, ccdrFit2graph);"
    "shd.val <- sapply(graph.path, pcalg::shd, g1);"
    "print(shd.val);"
    "z <- which(shd.val == min(shd.val));"
    "z <- z[length(z)];"
    "graph.shd <- permutenodes(graph.path[[z]], q);"
    "cg <- compare.graph(graph.shd, g);"
    "cg";
    
    Rcpp::NumericVector cg = R.parseEval(cmd);
    
    exit(0);
}


