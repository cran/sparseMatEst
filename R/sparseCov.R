
#' Sparse covariance matrix estimator with error control
#'
#'  Given a data matrix, \code{sparseCov} estimates the covariance
#'  matrix for the data under the assumption that the true covariance
#'  matrix is sparse, i.e. that most of the off-diagonal entries are
#'  equal to zero.
#'
#'  The algorithm begins with the empirical covariance estimator as
#'  computed by \code{cov}.  It then iteratively computes covariance
#'  estimators with false positive rates of \code{alf}, \code{alf}^2,
#'  and so on until \code{iter} estimators have been constructed.
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
#' @param alf  false positive rate in [0,1], Default is 0.5
#' @param iter number of iterates, Default is 10
#' @param pnrm norm to use, Default = Inf
#' @param THRSH Type of thresholding used;
#'              Takes values: hard, soft, adpt, scad.
#'
#' @return a list of arrays containing \code{iter+1} sparse
#'  covariance matrices corresponding to false positive rates of
#'  \code{1}, \code{alf},
#'  \code{alf}^2,..., \code{alf}^\code{iter}.  Each list
#'  corresponds to one type of thresholding chosen by \code{THRSH}.
#'
#' @author
#'   Adam B Kashlak \email{kashlak@ualberta.ca}
#'
#' @references
#'
#'   Kashlak, Adam B., and Linglong Kong.
#'   "A concentration inequality based methodology
#'   for sparse covariance estimation." arXiv
#'   preprint arXiv:1705.02679 (2017).
#'
#' @examples
#' # Generate four sparse covariance matrix estimators
#' # with false positive rates of 0.5, 0.25, 0.125,
#' # and 0.0625
#' n = 30
#' k = 50
#' dat = matrix(rnorm(n*k),n,k)
#' out = sparseCov( dat, alf=0.5, iter=4, THRSH=c("hard","soft") )
#'   par(mfcol=c(2,2))
#'   lab = c(1,0.5,0.5^2,0.5^3,0.5^4);
#'   for( i in 2:5 )
#'     image( out$hard[,,i]!=0, main=lab[i] )
#'   for( i in 2:5 )
#'     image( log(abs(out$hard[,,i] - out$soft[,,i])), main=lab[i] )
#' @export

sparseCov <- function(
  dat, alf=0.5, iter=10, pnrm=Inf, THRSH='hard'
){
  return(
    sparseMat(
      cov(dat),ncol(dat),alf,iter,pnrm,THRSH
    )
  );
}
