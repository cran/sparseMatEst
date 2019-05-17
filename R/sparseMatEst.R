#' Sparse Matrix Estimation and Inference
#'
#' The \code{sparseMatEst} library provides functions for estimating
#' sparse covariance and precision matrices with error control.
#'
#' Given a data matrix, this package contains two main functions
#' used to estimate the covariance matrix and the precision matrix
#' of the data under the assumption of sparsity, that most
#' off-diagonal entries are zero.  This is achieved by selecting
#' a false positive rate corresponding to the probability that
#' a true zero entry is falsely chosen to be non-zero by the estimator.
#'
#' The false positive rate can be treated as an interpretable
#' penalization parameter.  Setting this to zero will return
#' a diagonal matrix.  Choosing a false positive rate away from
#' zero will allow for the estimator to contain some non-zero
#' off-diagonal entries.
#'
#' Future updates coming Fall 2019 include inferential tools
#' based on these
#' sparse matrix estimators.  These include a variant of
#' linear and quadratic discriminant analysis, fitting a
#' Gaussian mixture assuming sparsity, network clustering
#' algorithm, and a method to fit a random design linear
#' regression model.
#'
#' @author
#' Adam B Kashlak \email{kashlak@ualberta.ca}
#'
#' @references
#'
#'   Kashlak, Adam B., and Linglong Kong.
#'   "A concentration inequality based methodology
#'   for sparse covariance estimation." arXiv
#'   preprint arXiv:1705.02679 (2017).
#'
#'   Kashlak, Adam B. "Non-asymptotic error controlled
#'   sparse high dimensional precision matrix estimation."
#'   arXiv preprint arXiv:1903.10988 (2019).
#'
#' @importFrom stats cov quantile rbinom runif
#' @importFrom glasso glasso
#'
"_PACKAGE"
#> [1] "_PACKAGE"
