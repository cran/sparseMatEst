#'  Sparse precision matrix estimator with error control
#'
#'  Given a data matrix, \code{sparsePrec} estimates the precision
#'  matrix for the data under the assumption that the true precision
#'  matrix is sparse, i.e. that most of the off-diagonal entries are
#'  equal to zero.
#'
#'  The algorithm begins with an initial precision matrix estimator
#'  being either \code{regMeth} = 'glasso' for debiased graphical lasso
#'  or \code{regMeth} = 'ridge' for debiased ridge estimator.
#'  It then iteratively computes covariance
#'  estimators with false positive rates of \code{alf}, \code{alf}^2,
#'  and so on until \code{iter} estimators have been constructed.
#'
#'  If \code{dat} is not included, but \code{prec} is, then the
#'  code runs as before but using the argument \code{prec} as the
#'  initial precision matrix estimator.  In this case, the \code{regMeth}
#'  is not considered.
#'
#'  The norm chosen determines the topology on the space of matrices.
#'  \code{pnorm} defaults to \code{Inf} being the operator norm
#'  or maximal eigenvalue.  This is theoretically justified to work
#'  in the references.
#'  Other norms could be considered, but their performance is not
#'  a strong.
#'
#'  Four thresholding methods are implemented.  \code{THRSH} defaults
#'  to hard thresholding where small matrix entries are set to zero
#'  while large entries are not affected.  The soft and adpt thresholds
#'  shrink all entries towards zero.  The scad threshold interpolates
#'  between hard and soft thresholding.  More details can be found
#'  in the references.
#'
#' @param dat  nxk data matrix, n observations, k dimensions
#' @param prec start with a precision estimator instead of \code{dat}
#' @param alf  false positive rate in [0,1], Default is 0.5
#' @param iter number of iterates, Default is 10
#' @param pnrm norm to use, Default = Inf
#' @param THRSH Type of thresholding used;
#'              Takes values: hard, soft, adpt, scad.
#' @param rho   Penalization parameter for the initial estimator, Default is 1
#' @param regMeth Type of initial estimator, Default is 'glasso'
#'
#' @return a list of arrays containing \code{iter+1} sparse
#'  precision matrices corresponding to false positive rates of
#'  \code{1}, \code{alf},
#'  \code{alf}^2,..., \code{alf}^\code{iter}.  Each list
#'  corresponds to one type of thresholding chosen by \code{THRSH}.
#'
#' @author
#'   Adam B Kashlak \email{kashlak@ualberta.ca}
#'
#' @references
#'
#'   Kashlak, Adam B. "Non-asymptotic error controlled
#'   sparse high dimensional precision matrix estimation."
#'   arXiv preprint arXiv:1903.10988 (2019).
#'
#' @examples
#' # Generate four sparse covariance matrix estimators
#' # with false positive rates of 0.5, 0.25, 0.125,
#' # and 0.0625
#' n = 30
#' k = 50
#' dat = matrix(rnorm(n*k),n,k)
#' out = sparsePrec( dat, alf=0.5, iter=4, THRSH=c("hard","soft") )
#'   par(mfcol=c(2,2))
#'   lab = c(1,0.5,0.5^2,0.5^3,0.5^4);
#'   for( i in 2:5 )
#'     image( out$hard[,,i]!=0, main=lab[i] )
#'   for( i in 2:5 )
#'     image( log(abs(out$hard[,,i] - out$soft[,,i])), main=lab[i] )
#' @export


sparsePrec <- function(
  dat=0,prec=0, alf=0.5, iter=10, pnrm=Inf, THRSH='hard',
  rho=1, regMeth='glasso'
){
  if(length(dat)>1){
    return(
      sparseMat(
        desparsePrec(dat,rho,type=regMeth),ncol(dat),alf,iter,pnrm,THRSH
      )
    );
  } else if(length(prec)>1){
    return(
      sparseMat(
        prec,ncol(prec),alf,iter,pnrm,THRSH
      )
    );
  } else {
    warning("ERROR: No data or precision matrix!");
  }
}
