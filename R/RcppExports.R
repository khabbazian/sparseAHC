#' agglomerative hierarchical clustering for sparse similarity matrices
#'
#'@param S sparse similarity matrix 
#'@param linkage is c(average, single, complete) 
#'@param check.symmetry logical. If TRUE, checks if the input similarity matrix is symmetric. 
#'@param noOrder logical. If TRUE, does not return the order of points which is required for plotting dendrogram to save memory and time
#'
#'@return hclust object
#'
#'@examples
#'library("Matrix");
#'library("igraph");
#'library("sparseHC");
#'
#'A = Matrix(0, nrow = 6, ncol = 6, sparse = TRUE);
#'A[1,2]= 2; A[1,5]= 3; A[1,4]= 2.5; A[2,5]= 4;
#'A[3,4]= 2; A[3,5]= 6; A[4,5]= 5; A[5,6]= 4.6;
#'A = A + t(A);
#'
#'G = graph.adjacency(A, mode = "undirected", weighted=TRUE);
#'plot(G,edge.label=E(G)$weight,vertex.label=V(G)-1)
#'
#'H=sparseHC(A,"average",FALSE)
#'plot(H);
#' 
#'@export
sparseHC <- function(S, linkage=c("average", "single", "complete"), check.symmetry=TRUE, noOrder = FALSE) {
    
    linkage = match.arg(linkage);
    if(!isSymmetric(S)) stop("the input similarity matrix has to be symmetric")
    if( !is(S, "dgCMatrix") ){
        warning(paste0("Input matrix type changed from ", class(S), " to dgCMatrix"));
        S  <-  as(S, "dgCMatrix");
    }
    stopifnot( is(S, "dgCMatrix") )
    .Call('sparseHC_run_sparseHC', PACKAGE = 'sparseHC', S, linkage, noOrder)
}

