#' Steihaug truncated-CG trust-region Newton minimization
#'
#' Thread-safe C++ port of the trust-region Newton minimizer with the
#' Steihaug truncated conjugate-gradient subproblem from the Rust crate
#' \href{https://github.com/jolars/basin}{basin} (\code{TrustRegion} with
#' its default \code{Steihaug} subproblem, driven by basin's executor). Each
#' iteration approximately minimizes the quadratic model
#' \eqn{m(p) = f + g^T p + \frac{1}{2} p^T B p}{m(p) = f + g'p + p'Bp/2}
#' over \eqn{\|p\| \le \Delta}{||p|| <= Delta} by conjugate gradients,
#' stopping at the boundary on negative curvature. Unlike [trust()] it never
#' factorizes the Hessian, and with \code{hessvec} it never forms it at all
#' (basin's matrix-free mode).
#'
#' Evaluation pattern (identical to basin): \code{fn} and \code{gr} once at
#' the start; \code{hess} once per outer iteration (exact mode) or
#' \code{hessvec} once per CG product (matrix-free mode); \code{fn} at each
#' trial point; \code{gr} after each accepted step. \code{fn} may return
#' \code{Inf} (or \code{NaN}) at a trial point to reject it.
#'
#' @param par numeric vector of starting parameters.
#' @param fn objective function of the parameter vector returning a scalar.
#' @param gr gradient function returning a numeric vector of length
#'   \code{length(par)}.
#' @param hess Hessian function returning a symmetric
#'   \code{length(par)} by \code{length(par)} matrix. May be \code{NULL} when
#'   \code{hessvec} is supplied.
#' @param hessvec optional Hessian-vector product \code{function(x, v)}
#'   returning \eqn{H(x) v}{H(x) \%*\% v}. When supplied, the matrix-free mode
#'   is used and \code{hess} is never called.
#' @param control a list of control settings (defaults in brackets; a
#'   \code{NULL} tolerance or budget is disabled):
#'   \describe{
#'     \item{rinit}{initial trust region radius [1]}
#'     \item{rmax}{maximum trust region radius [100]}
#'     \item{eta}{step acceptance threshold on the reduction ratio, in
#'       \eqn{[0, 1/4)} [0.125]}
#'     \item{maxInner}{radius reductions (re-solves) per outer iteration [10]}
#'     \item{cgMaxit}{CG iteration cap per subproblem; 0 means
#'       \code{length(par)} [0]}
#'     \item{kappa, theta}{CG forcing rule: stop when
#'       \eqn{\|r\| < \min(\kappa, \|g\|^\theta)\|g\|}{||r|| < min(kappa, ||g||^theta) ||g||}
#'       [0.5, 0.5]}
#'     \item{maxit}{maximum outer iterations [100]}
#'     \item{maxCostEvals, maxGradEvals}{evaluation budgets [NULL]}
#'     \item{targetCost}{stop once the objective is at or below this [NULL]}
#'     \item{gradTol}{absolute gradient tolerance on
#'       \eqn{\|g\|_2}{||g||_2} [1e-8]}
#'     \item{relGradTol}{gradient tolerance relative to the initial gradient
#'       norm [NULL]}
#'     \item{stepTol, relStepTol}{absolute / relative step tolerance [NULL]}
#'     \item{costTol, relCostTol}{absolute / relative objective change
#'       tolerance [NULL]}
#'     \item{trace}{if \code{TRUE}, return a per-attempt trace [FALSE]}
#'   }
#' @param ... further arguments passed to \code{fn}, \code{gr}, \code{hess}
#'   and \code{hessvec}.
#'
#' @return A list with components \code{par}, \code{value}, \code{gradient},
#'   \code{iterations}, \code{converged}, \code{status} (integer),
#'   \code{message} (the basin termination reason, e.g.
#'   \code{"GradientTolerance"}, \code{"MaxIter"}), \code{radius} (final
#'   trust radius), \code{counts} (evaluations of \code{fn}, \code{gr},
#'   \code{hess}, \code{hessvec}), and, when \code{control$trace} is
#'   \code{TRUE}, \code{trace}, a data frame with one row per subproblem
#'   attempt.
#' @references Steihaug, T. (1983). The conjugate gradient method and trust
#'   regions in large scale optimization. SIAM Journal on Numerical Analysis,
#'   20(3), 626-637.
#'
#'   Nocedal, J. and Wright, S. J. (2006). Numerical Optimization (2nd ed.),
#'   Algorithms 4.1 and 7.2. Springer.
#' @seealso [trust()], [bobyqa()], [newuoa()]
#' @examples
#' fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
#' grr <- function(x) c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
#'                      200 * (x[2] - x[1]^2))
#' hr <- function(x) matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
#'                            -400 * x[1], 200), 2, 2)
#' steihaug(c(-1.2, 1), fr, grr, hr)
#'
#' # matrix-free: only Hessian-vector products
#' hv <- function(x, v) drop(hr(x) %*% v)
#' steihaug(c(-1.2, 1), fr, grr, hessvec = hv)$counts
#' @export
steihaug <- function(par, fn, gr, hess = NULL, hessvec = NULL,
                     control = list(), ...) {
  if (!is.numeric(par) || length(par) < 1) stop("par must be a non-empty numeric vector")
  if (!all(is.finite(par))) stop("par not all finite")
  if (!is.function(fn)) stop("fn must be a function")
  if (!is.function(gr)) stop("gr must be a function")
  if (is.null(hess) && is.null(hessvec))
    stop("one of 'hess' or 'hessvec' must be supplied")
  if (!is.null(hess) && !is.function(hess)) stop("hess must be a function")
  if (!is.null(hessvec) && !is.function(hessvec))
    stop("hessvec must be a function")

  ctrl <- list(rinit = 1, rmax = 100, eta = 0.125, maxInner = 10L,
               cgMaxit = 0L, kappa = 0.5, theta = 0.5, maxit = 100L,
               maxCostEvals = NULL, maxGradEvals = NULL, targetCost = NULL,
               gradTol = 1e-8, relGradTol = NULL, stepTol = NULL,
               relStepTol = NULL, costTol = NULL, relCostTol = NULL,
               trace = FALSE)
  if (!is.list(control)) stop("control must be a list")
  unused <- setdiff(names(control), names(ctrl))
  if (length(unused))
    warning("unused control arguments ignored: ",
            paste(unused, collapse = ", "))
  for (nm in intersect(names(control), names(ctrl))) {
    ctrl[nm] <- list(control[[nm]])
  }

  for (nm in c("rinit", "rmax", "eta", "kappa", "theta")) {
    if (!is.numeric(ctrl[[nm]]) || length(ctrl[[nm]]) != 1 ||
        !is.finite(ctrl[[nm]]))
      stop("control$", nm, " must be a single finite number")
  }
  for (nm in c("maxInner", "cgMaxit", "maxit", "maxCostEvals", "maxGradEvals")) {
    v <- ctrl[[nm]]
    if (is.null(v) && nm %in% c("maxCostEvals", "maxGradEvals")) next
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 0 ||
        v > .Machine$integer.max)
      stop("control$", nm, " must be a single non-negative finite count")
  }
  for (nm in c("targetCost", "gradTol", "relGradTol", "stepTol", "relStepTol",
               "costTol", "relCostTol")) {
    v <- ctrl[[nm]]
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1 || !is.finite(v)))
      stop("control$", nm, " must be NULL or a single finite number")
  }

  nn <- names(par)
  dots <- list(...)
  wrap1 <- function(f) {
    force(f)
    function(x) {
      names(x) <- nn
      do.call(f, c(list(x), dots))
    }
  }
  hv1 <- NULL
  if (!is.null(hessvec)) {
    hv1 <- function(x, v) {
      names(x) <- nn
      do.call(hessvec, c(list(x, v), dots))
    }
  }
  ret <- .Call(`_RcppTrust_steihaug_solve_r`, as.double(par), wrap1(fn),
               wrap1(gr), if (is.null(hess)) NULL else wrap1(hess), hv1, ctrl)
  names(ret$par) <- nn
  names(ret$gradient) <- nn
  ret
}
