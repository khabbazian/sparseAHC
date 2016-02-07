#' agglomerative hierarchical clustering for sparse similarity matrices
#'
#'@param S sparse similarity matrix of type dgCMatrix
#'@param linkage is c(average, single, complete) 
#'@param check.symmetry logical. If TRUE, checks if the input similarity matrix is symmetric. 
#'
#'@return hclust object
#'
#'@examples
#'library("Matrix")
#'library("igraph")
#'
#'A <- Matrix(0, nrow = 6, ncol = 6, sparse = TRUE)
#'A[1,2] <- 2; A[1,5] <- 3; A[1,4] <- 2.5; A[2,5] <- 4
#'A[3,4] <- 2; A[3,5] <- 6; A[4,5] <- 5; A[5,6] <- 4.6
#'A <- A + t(A)
#'
#'G <- graph.adjacency(A, mode = "undirected", weighted=TRUE)
#'plot(G,edge.label=E(G)$weight,vertex.label=V(G)-1)
#'
#'H <- sparseAHC(A, "average", TRUE)
#'plot(H)
#' 
#'@export
sparseAHC <- function(S, linkage=c("average", "single", "complete"), check.symmetry=TRUE) {
    
    linkage = match.arg(linkage)

    if( !is(S, "dgCMatrix") ){
        warning(paste0("Input matrix type changed from ", class(S), " to dgCMatrix"))
        S  <- as(S, "dgCMatrix")
    }
    if( check.symmetry ){
        if ( ! .Call('sparseAHC_dgCIsSymmetric', PACKAGE = 'sparseAHC', S, 100*.Machine$double.eps) ){
            stop("The input similarity matrix should be symmetric.")
        }
    }

    stopifnot( is(S, "dgCMatrix") )
    .Call('sparseAHC_run_sparseAHC', PACKAGE = 'sparseAHC', S, linkage)
}

#'Test if a Matrix of type dgCMatrix is Symmetric
#'Test if a Matrix of type dgCMatrix is Symmetric
#'
#'@param object any R object; a ‘matrix’ for the matrix method.
#'@param tol numeric scalar >= 0.  Smaller differences are not considered.
#'@param ... further arguments passed to methods; the matrix method passes.
#'
#'@export
isSymmetric.dgCMatrix <- function(object, tol=100*.Machine$double.eps, ...){
    .Call('sparseAHC_dgCIsSymmetric', PACKAGE = 'sparseAHC', S = object, eps = tol)
}

