#
#  ccdr-utils.R
#  ccdri
#
#  Created by Bryon Aragam (local) on 5/28/14.
#  Copyright (c) 2014-2015 Bryon Aragam (local). All rights reserved.
#

#
# PACKAGE CCDR: Utils
#
#   CONTENTS:
#     check_if_matrix
#     list_classes
#     check_list_class
#     col_classes
#     cor_vector
#     cor_vector_intervention
#

# Special function to check if an object is EITHER matrix or Matrix object
check_if_matrix <- function(m){
    is.matrix(m) || inherits(m, "Matrix")
} # END .CHECK_IF_MATRIX

check_if_data_matrix <- function(df){
    is.data.frame(df) || is.matrix(df)
} # END .CHECK_IF_DATA_MATRIX

# Count missing values in a matrix or data.frame
count_nas <- function(df){
    if( !check_if_data_matrix(df)){
        stop("Input must be a data.frame or a matrix!")
    }

    sum(is.na(df))
} # END .COUNT_NAS

# Special function to return types for each element in a list
list_classes <- function(li){
    unlist(lapply(li, class))
} # END .LIST_CLASSES

# Return TRUE if every element of li inherits check.class, FALSE otherwise
check_list_class <- function(li, check.class){
    if(length(li) == 0){
        warning("List contains no elements!")

        TRUE # default to true if empty
    }

    all(unlist(lapply(li, function(x) inherits(x, check.class))))
} # END .CHECK_LIST_CLASS

# Output the class of each column in X, return as a character vector
col_classes <- function(X){
    if( !check_if_data_matrix(X)){
        stop("Input must be a data.frame or a matrix!")
    }

    apply(X, 2, class)
} # END .COL_CLASSES

cor_vector <- function(X){
# This is now implicitly checked via col_classes
#     if( !is.data.frame(X) && !is.matrix(X)){
#         stop("Input must either be a data.frame or a matrix!")
#     }

    check.numeric <- (col_classes(X) != "numeric")
    if( any(check.numeric)){
        not.numeric <- which(check.numeric)
        stop(paste0("Input columns must be numeric! Columns ", paste(not.numeric, collapse = ", "), " are non-numeric."))
    }

    if( any(dim(X) < 2)){
        stop("Input must have at least 2 rows and columns!") # 2-8-15: Why do we check this here?
    }

    cors <- cor(X)
    cors <- cors[upper.tri(cors, diag = TRUE)]

    cors
} # END .COR_VECTOR

cor_vector_intervention <- function(X, intervention = NULL) {
## Similar to cor_vector
## Now there are interventions
    check.numeric <- (col_classes(X) != "numeric")
    if( any(check.numeric)){
        not.numeric <- which(check.numeric)
        stop(paste0("Input columns must be numeric! Columns ", paste(not.numeric, collapse = ", "), " are non-numeric."))
    }

    if( any(dim(X) < 2)){
        stop("Input must have at least 2 rows and columns!") # 2-8-15: Why do we check this here?
    }

    pp <- ncol(X)
    if(!is.null(intervention)) {
        if(length(intervention) != nrow(X)) stop("Intervention size does not match")
    } else intervention <- as.integer(rep(pp + 1, nrow(X)))

    ivj <- sort(unique(intervention)) # all the j's that has intervention (including pp+1)
    len <- length(ivj)
    cors <- vector("list", len)

    indexj <- as.integer(rep(len - 1, pp + 1))
    if(len > 1) for(j in 1:(len - 1)) {
        jj <- ivj[j]
        indexj[jj] <- j - 1
        corsjj <- cor(X[intervention != jj, ])
        cors[[j]] <- corsjj[upper.tri(corsjj, diag = TRUE)]
    }
    corsjj <- cor(X[intervention == pp + 1, ])
    cors[[len]] <- corsjj[upper.tri(corsjj, diag = TRUE)]
    cors <- unlist(cors)
    return(list(cors = cors, indexj = indexj))
} # END .COR_VECTOR_INTERVETION
