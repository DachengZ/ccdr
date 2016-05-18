### Ensure get.rhos function is available to all tests
get.rhos <- function(x) UseMethod("get.rhos", x)

get.rhos.ccdrFit <- function(cf){
    cf$sbm$sigmas
}
