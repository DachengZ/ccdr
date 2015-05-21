<!-- README.md is generated from README.Rmd. Please edit that file -->
ccdr
====

`ccdr` implements the CCDr structure learning algorithm described in [Aragam and Zhou (2015), Journal of Machine Learning Research](http://arxiv.org/abs/1401.0852). Based on observational data, this algorithm estimates the structure a Bayesian network (aka edges in a DAG) using penalized maximum likelihood based on L1 or concave (MCP) regularization.

Presently, this package consists of a single method that implements the main algorithm; more functionality will be provided in the future. To generate data from a given Bayesian network and/or simulate random networks, the following R packages are recommended:

-   `bnlearn`: [www.bnlearn.com](http://www.bnlearn.com), [bnlearn on CRAN](http://cran.r-project.org/web/packages/bnlearn/index.html)
-   `pcalg`: [pcalg on CRAN](http://cran.r-project.org/web/packages/pcalg/index.html)
-   `igraph`: [<http://igraph.org/r/>](http://igraph.org/r/), [igraph on CRAN](http://cran.r-project.org/web/packages/igraph/index.html)

Installation
============

Please note that this package is currently in a development state. In order to install the package, you will need both `devtools` and `Rcpp`. We are still working out the kinks of getting the package to compile on different systems, so please [let us know](https://github.com/itsrainingdata/ccdr/issues) if you have any issues.

If you have never installed packages using `Rcpp` before, we recommend checking out the following resources which should get you started:

-   [Rcpp for Seamless R and C++ Integration](http://www.rcpp.org)
-   [Rcpp: Seamless R and C++ Integration](http://www.jstatsoft.org/v40/i08/), Dirk Eddelbuettel, Romain Francois, Journal of Statistical Software, Vol. 40, Issue 8, Apr 2011.
-   [High performance functions with Rcpp](http://adv-r.had.co.nz/Rcpp.html)

The code is hosted on github and can be installed using `devtools`:

``` r
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
```

**NOTE:** Windows users will need to make sure that Rtools is installed in order to build this package. To check if you have Rtools installed, you can use `devtools::find_rtools()`. For more details: [<http://cran.r-project.org/bin/windows/Rtools/>](http://cran.r-project.org/bin/windows/Rtools/).

If you find any bugs, please report them here on [github](https://github.com/itsrainingdata/ccdr/issues).

Examples
========

The basic usage is as follows:

``` r
### Specify dimensions
nn <- 20
pp <- 100

### Generate a random Gaussian matrix
dat <- matrix(rnorm(nn * pp), nrow = nn)

### Run the ccdr algorithm
ccdr.run(data = dat)
```

This trivial example uses uncorrelated normal data, which is not very interesting. In order to do some interesting calculations, we first need to generate data according to some pre-specified DAG structure. The `ccdr` package does not provide this functionality: You can use either the `bnlearn` package or the `pcalg` package to generate random DAGs and random data.

Example using `bnlearn`
-----------------------

Example using `pcalg`
---------------------

The `pcalg` package provides two useful functions: `randomDAG` for generating a random *ordered* DAG and `rmvDAG` for generating random data according to the structural equation model implied by the DAG.

**NOTE:** Since `randomDAG` produces an ordered DAG, we should permute the columns in the simulated dataset in order to obfuscate this information. If we know this ordering in advance, there are much better methods for estimating the DAG structure using autoregressive models.

``` r
library("pcalg")

### Set up the model parameters
nn <- 20                # How many samples to draw?
pp <- 100               # How many nodes in the DAG?
num.edges <- 100        # How many *expected* edges in the DAG?
ss <- num.edges / pp    # This is the expected number of parents *per node*

### Generate a random DAG using the pcalg method randomDAG
beta.min <- 0.5
beta.max <- 2
edge.pr <- 2 * ss / (pp - 1)
g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!

### Generate random data according to this DAG using the method rmvDAG
dat <- rmvDAG(n = nn, dag = g, errDist = "normal")
pi <- sample(1:pp)
dat <- dat[, pi] # permute the columns to randomize node ordering

### Run the algorithm
final2 <- ccdr.run(data = dat, lambdas.length = 10, alpha = 10, verbose = FALSE)
```