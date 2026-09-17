# R shims for the thread-safe C++ ports of Powell's BOBYQA and NEWUOA.
# These reproduce minqa 1.2.8's R/minqa.R argument handling,
# defaults, warnings and return values line for line; only the .Call
# target differs (the ported solvers in src/minqa_*.cpp instead of the
# Fortran).

# minqa's commonArgs(): establishes control defaults, assigning `n` and an
# attribute-free `par` into the caller's frame exactly as minqa does.
.minqaCommonArgs <- function(par, fn, ctrl, rho) {
  rho$n <- n <- length(rho$par <- as.double(par))
  stopifnot(all(is.finite(par)),
            is.function(fn),
            length(formals(fn)) >= 1)
  cc <- do.call(function(npt = min(n+2L, 2L * n), rhobeg = NA,
                         rhoend = NA, iprint = 0L, maxfun=10000L,
                         obstop=TRUE, force.start=FALSE,...) {
    if (length(list(...))>0) warning("unused control arguments ignored")
    list(npt = npt, rhobeg = rhobeg, rhoend = rhoend,
         iprint = iprint, maxfun = maxfun, obstop = obstop,
         force.start = force.start)
  }, ctrl)

  ctrl <- new.env(parent = emptyenv())
  lapply(names(cc), function(nm) assign(nm, cc[[nm]], envir = ctrl))

  ctrl$npt <- as.integer(max(n + 2L, min(ctrl$npt, ((n+1L)*(n+2L)) %/% 2L)))
  if (ctrl$npt > (2 * n + 1))
    warning("Setting npt > 2 * length(par) + 1 is not recommended.")

  if (is.na(ctrl$rhobeg))
    ctrl$rhobeg <- min(0.95, 0.2 * max(abs(par)))
  if (is.na(ctrl$rhoend)) ctrl$rhoend <- 1.0e-6 * ctrl$rhobeg
  stopifnot(0 < ctrl$rhoend, ctrl$rhoend <= ctrl$rhobeg)

  if (ctrl$maxfun < 10 * n^2)
    warning("maxfun < 10 * length(par)^2 is not recommended.")
  ctrl
}

#' Nonlinear optimization with box constraints (BOBYQA)
#'
#' Thread-safe C++ port of \code{minqa::bobyqa()}: minimizes a function of
#' many variables subject to box constraints by a trust region method that
#' forms quadratic models by interpolation, using M. J. D. Powell's BOBYQA
#' algorithm. Arguments, defaults, and return value are those of
#' \code{minqa::bobyqa()}; see its documentation for details.
#'
#' @param par numeric vector of starting parameters.
#' @param fn function to be minimized; its first argument must be the
#'   parameter vector and it must return a scalar numeric value.
#' @param lower,upper numeric vectors of lower and upper bounds (recycled
#'   when of length 1).
#' @param control a list of control settings: \code{npt}, \code{rhobeg},
#'   \code{rhoend}, \code{iprint}, \code{maxfun}, \code{obstop} and
#'   \code{force.start}, as in \code{minqa::bobyqa()}.
#' @param ... further arguments passed to \code{fn}.
#'
#' @return A list of class \code{c("bobyqa", "minqa")} with components
#'   \code{par}, \code{fval}, \code{feval}, \code{ierr} and \code{msg}.
#' @references M. J. D. Powell (2009), "The BOBYQA algorithm for bound
#'   constrained optimization without derivatives", Report No. DAMTP
#'   2009/NA06, Centre for Mathematical Sciences, University of Cambridge.
#' @seealso [newuoa()], and [minqa_c_api] for calling the
#'   solvers from parallel C/C++ code.
#' @examples
#' fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
#' bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4))
#' @export
bobyqa <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
{
  nn <- names(par)
  ctrl <- .minqaCommonArgs(par, fn, control, environment())
  n <- length(par)
  fn1 <- function(x) {
    names(x) <- nn
    fn(x, ...)
  }
  checkObj <- fn1(par)
  if(length(checkObj) > 1 || !is.numeric(checkObj))
    stop("Objective function must return a single numeric value.")
  lower <- as.double(lower); upper <- as.double(upper)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  stopifnot(length(lower) == n, length(upper) == n, all(lower < upper))
  if (any(par < lower | par > upper)) {
    if (ctrl$obstop)
      stop("Starting values violate bounds")
    else {
      # kept verbatim from minqa (including its pmax(par, upper))
      par <- pmax(lower, pmax(par, upper))
      warning("Some parameters adjusted to nearest bound")
    }
  }
  rng <- upper - lower

  if (any(rng < 2 * ctrl$rhobeg)) {
    warning("All upper - lower must be >= 2*rhobeg. Changing rhobeg")
    ctrl$rhobeg <- 0.2 * min(rng)
  }

  verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
  if (all(is.finite(upper)) && all(is.finite(lower)) &&
      all(par >= lower) && all(par <= upper) ) {
    if (verb) cat("ctrl$force.start = ", ctrl$force.start,"\n")
    if (!ctrl$force.start) {
      i <- rng < ctrl$rhobeg
      if (any(i)) {
        par[i] <- lower[i] + ctrl$rhobeg
        warning("Some parameters adjusted away from lower bound")
      }
      i <- rng < ctrl$rhobeg
      if (any(i)) {
        par[i] <- upper[i] - ctrl$rhobeg
        warning("Some parameters adjusted away from upper bound")
      }
    }
  }
  if (verb) {
    cat("npt =", ctrl$npt, ", n = ",n,"\n")
    cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
  }
  if(ctrl$iprint > 0)
    cat("start par. = ", par, "fn = ", checkObj, "\n")

  retlst <- .Call(`_RcppTrust_minqa_bobyqa_r`, as.double(par), lower, upper,
                  ctrl, fn1)
  if (retlst$ierr > 0){
    if (retlst$ierr == 10) {
      retlst$ierr<-2
      retlst$msg<-"bobyqa -- NPT is not in the required interval"
    } else if (retlst$ierr == 320) {
      retlst$ierr<-5
      retlst$msg<-"bobyqa detected too much cancellation in denominator"
    } else if (retlst$ierr == 390) {
      retlst$ierr<-1
      retlst$msg<-"bobyqa -- maximum number of function evaluations exceeded"
    } else if (retlst$ierr == 430) {
      retlst$ierr<-3
      retlst$msg<-"bobyqa -- a trust region step failed to reduce q"
    } else if (retlst$ierr == 20) {
      retlst$ierr<-4
      retlst$msg<-"bobyqa -- one of the box constraint ranges is too small (< 2*RHOBEG)"
    }
  } else {
    retlst$msg<-"Normal exit from bobyqa"
  }
  retlst
}

#' Unconstrained derivative-free optimization (NEWUOA)
#'
#' Thread-safe C++ port of \code{minqa::newuoa()}: minimizes a function of
#' many variables by a trust region method that forms quadratic models by
#' interpolation, using M. J. D. Powell's NEWUOA algorithm.
#'
#' @inheritParams bobyqa
#' @param control a list of control settings: \code{npt}, \code{rhobeg},
#'   \code{rhoend}, \code{iprint} and \code{maxfun}, as in
#'   \code{minqa::newuoa()}.
#'
#' @return A list of class \code{c("newuoa", "minqa")} with components
#'   \code{par}, \code{fval}, \code{feval}, \code{ierr} and \code{msg}.
#' @references M. J. D. Powell (2006), "The NEWUOA software for
#'   unconstrained optimization without derivatives", in Large-Scale
#'   Nonlinear Optimization, Springer, 255-297.
#' @seealso [bobyqa()], [minqa_c_api]
#' @examples
#' fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
#' newuoa(c(1, 2), fr)
#' @export
newuoa <- function(par, fn, control = list(), ...)
{
  nn <- names(par)
  ctrl <- .minqaCommonArgs(par + 0, fn, control, environment())
  n <- length(par)
  fn1 <- function(x) {
    names(x) <- nn
    fn(x, ...)
  }
  checkObj <- fn1(par)
  if(length(checkObj) > 1 || !is.numeric(checkObj))
    stop("Objective function must return a single numeric value.")
  verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
  if (verb) {
    cat("npt =", ctrl$npt, ", n = ",n,"\n")
    cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
  }
  if(ctrl$iprint > 0)
    cat("start par. = ", par, "fn = ", checkObj, "\n")

  retlst <- .Call(`_RcppTrust_minqa_newuoa_r`, as.double(par), ctrl, fn1)
  if (retlst$ierr > 0){
    if (retlst$ierr == 10) {
      retlst$ierr<-2
      retlst$msg<-"newuoa -- NPT is not in the required interval"
    } else if (retlst$ierr == 320) {
      retlst$ierr<-5
      retlst$msg<-"newuoa detected too much cancellation in denominator"
    } else if (retlst$ierr == 390) {
      retlst$ierr<-1
      retlst$msg<-"newuoa -- maximum number of function evaluations exceeded"
    } else if (retlst$ierr == 3701) {
      retlst$ierr<-3
      retlst$msg<-"newuoa -- a trust region step failed to reduce q"
    }
  } else {
    retlst$msg<-"Normal exit from newuoa"
  }
  retlst
}

#' Print method for minqa objects
#'
#' @param x an object inheriting from class \code{minqa}.
#' @param digits unused; kept for compatibility with \code{minqa}.
#' @param ... unused.
#' @return \code{invisible(x)}
#' @export
print.minqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("parameter estimates:", toString(x$par), "\n")
  cat("objective:", toString(x$fval), "\n")
  cat("number of function evaluations:", toString(x$feval), "\n")
  invisible(x)
}
