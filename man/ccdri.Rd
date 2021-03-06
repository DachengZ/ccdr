% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccdr.R
\docType{package}
\name{ccdri}
\alias{ccdri}
\alias{ccdri-package}
\title{ccdri: Structure learning for Bayesian networks using the CCDr algorithm.}
\description{
\code{ccdri} implements the CCDr structure learning algorithm for Bayesian networks.
}
\details{
Based on observational and interventional data, this algorithm estimates the structure of a Bayesian network
(aka edges in a DAG) using penalized maximum likelihood based on L1 or concave (MCP) regularization.
The main methods are:
\describe{
 \item{\code{\link{ccdr.run}}}{Executes the main algorithm.}
 \item{\code{\link{generate.lambdas}}}{Convenience method for generating a good sequence of regularization parameters.}
}
}

