
#' Sparse matrix generator
#'
#'  A way to generate sparse matrices for simulation and testing.
#'
#'  \code{genSparseMatrix} constructs a sparse matrix to be used for
#'  testing and simulation.  The \code{type} argument determines what
#'  type of matrix is produced while \code{k} effects the dimension and
#'  \code{rho} is an additional parameter to effect the output.
#'
#'  For \code{type = 'tri'}, a (k x k) tridiagonal matrix is returned
#'  with off-diagonal entries equal to \code{rho}.
#'
#'  For \code{type = 'arm'}, a (k x k) autoregressive matrix is returned
#'  with off-diagonal entries equal to \code{rho^|i-j|}.
#'
#'  For \code{type = 'band'}, a (k x k) banded matrix is returned with
#'  \code{rho} bands.
#'
#'  For \code{type = 'rand'}, a (k x k) matrix is returned with
#'  \code{rho} in (0,1) randomly selected off-diagonal entries set to
#'  be non-zero.
#'
#'  For \code{type = 'tree'}, the adjacency matrix for a k-deep binary tree
#'  is returned with off-diagonal entries set to \code{rho}.  The dimension
#'  of this matrix is k(k+1)/2.
#'
#'  For \code{type = 'multi'}, a (k x k) matrix is returned with \code{rho}
#'  off-diagonals set to ones.
#'
#'  For \code{type = 'block'}, a block diagonal matrix is returned with
#'  k blocks of size \code{rho}.
#'
#' @param type type of matrix to generate
#' @param k    dimension parameter
#' @param rho  additional parameter
#'
#' @return a matrix corresponding to the arguments chosen.
#'
#' @author
#'   Adam B Kashlak \email{kashlak@ualberta.ca}
#'
#' @examples
#'
#'   out = list();
#'   out[[1]]= genSparseMatrix( 20,0.5,"tri" );
#'   out[[2]]= genSparseMatrix( 20,0.5,"arm" );
#'   out[[3]]= genSparseMatrix( 20,5,"band" );
#'   out[[4]]= genSparseMatrix( 20,0.5,"rand" );
#'   out[[5]]= genSparseMatrix( 7,0.5,"tree" );
#'   out[[6]]= genSparseMatrix( 20,5,"multi" );
#'   out[[7]]= genSparseMatrix( 5,4,"block" );
#'
#'     par(mfrow=c(2,3));
#'     lab = c("tri","arm","band","rand","tree","multi","block");
#'     for( i in 2:7  )
#'       image(out[[i]],main=lab[i]);
#'
#' @export

genSparseMatrix <- function( k, rho, type ){
  switch( type,
    tri  = genTriDiag(k,rho),
    arm  = genArMat(k,rho),
    band = genBandedMat(k,rho),
    rand = genRandSpMat(k,rho),
    tree = genBinaryTreeMat(k,rho),
    multi= genMultDiag(k,1,rho),
    block= genBlockDiag(k,rho,1)
  )
}


