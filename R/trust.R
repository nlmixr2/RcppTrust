#' @useDynLib RcppTrust, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Internal helper called (via Rcpp::Function) once per objfun evaluation
# on the R-objfun path. Kept in R -- rather than reimplemented in C++ --
# specifically so that an error raised inside the user's objfun produces
# a genuine R "try-error" object (with its "condition" attribute), byte
# for byte what upstream trust()'s own try(objfun(theta, ...)) produces.
.trustTryEval <- function(f, theta) try(f(theta), silent = TRUE)

#' Positionally-indexed function pointer table for downstream packages
#'
#' Returns the list of native function pointers a consuming package (e.g.
#' nlmixr2est) resolves via the header-only registration pattern in
#' \code{inst/include/RcppTrust.h} -- the same pattern rxode2/n1qn1c/
#' lbfgsb3c use, so no linker dependency on this package's shared library
#' is required. Not intended to be called directly by end users.
#'
#' @return An unclassed list of external pointers.
#' @export
.RcppTrustPtr <- function() {
  .Call(`_RcppTrust_trust_ptr`, PACKAGE = "RcppTrust")
}

#' Non-Linear Optimization by a Trust Region Algorithm
#'
#' Thread-safe C++ port of \code{trust::trust()} (Geyer). Carries out a
#' minimization or maximization of a function using a trust region
#' algorithm; see the original package's documentation for the algorithm
#' and argument semantics, which this function reproduces exactly.
#'
#' @param objfun an R function computing value, gradient, and Hessian; see
#'   \code{trust::trust} for the exact contract.
#' @param parinit starting parameter values.
#' @param rinit starting trust region radius.
#' @param rmax maximum allowed trust region radius.
#' @param parscale optional parameter scale vector.
#' @param iterlim maximum number of iterations.
#' @param fterm function-value termination tolerance.
#' @param mterm model-decrease termination tolerance.
#' @param minimize if \code{TRUE} minimize, if \code{FALSE} maximize.
#' @param blather if \code{TRUE} return extra per-iteration information.
#' @param ... additional arguments passed to \code{objfun}.
#'
#' @return A list with the same components as \code{trust::trust()}.
#' @export
trust <- function(objfun, parinit, rinit, rmax, parscale,
    iterlim = 100, fterm = sqrt(.Machine$double.eps),
    mterm = sqrt(.Machine$double.eps),
    minimize = TRUE, blather = FALSE, ...) {
  if (!is.numeric(parinit)) stop("parinit not numeric")
  if (!all(is.finite(parinit))) stop("parinit not all finite")
  d <- length(parinit)

  if (missing(parscale)) {
    parscaleArg <- NULL
  } else {
    if (length(parscale) != d) stop("parscale and parinit not same length")
    if (!all(parscale > 0)) stop("parscale not all positive")
    if (!all(is.finite(parscale) & is.finite(1 / parscale)))
      stop("parscale or 1 / parscale not all finite")
    parscaleArg <- as.double(parscale)
  }
  if (!is.logical(minimize)) stop("minimize not logical")

  force(objfun)
  dotArgs <- list(...)
  objfun1 <- if (length(dotArgs)) {
    function(theta) do.call(objfun, c(list(theta), dotArgs))
  } else {
    function(theta) objfun(theta)
  }

  .Call(`_RcppTrust_trust_solve_r`, objfun1, .trustTryEval, as.double(parinit),
        as.double(rinit), as.double(rmax), parscaleArg, as.integer(iterlim),
        as.double(fterm), as.double(mterm), as.logical(minimize),
        as.logical(blather))
}
