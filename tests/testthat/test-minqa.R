# The BOBYQA/NEWUOA ports must reproduce minqa 1.2.8's results exactly
# (same evaluations, same iterates) for runs that never enter BOBYQA's
# RESCUE step, where minqa's rescue.f evaluates at uninitialized memory.

fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
extRosen <- function(x) {
  n <- length(x)
  sum(100 * (x[-1] - x[-n]^2)^2 + (1 - x[-n])^2)
}
chebyquad <- function(x) {
  n <- length(x)
  res <- numeric(n)
  for (i in 1:n) {
    rr <- 0
    for (k in 1:n) {
      z7 <- 1; z2 <- 2 * x[k] - 1; z8 <- z2; j <- 1
      while (j < i) { z6 <- z7; z7 <- z8; z8 <- 2 * z2 * z7 - z6; j <- j + 1 }
      rr <- rr + z8
    }
    rr <- rr / n
    if (2 * trunc(i / 2) == i) rr <- rr + 1 / (i * i - 1)
    res[i] <- rr
  }
  sum(res * res)
}

# Bitwise agreement is established on x86_64, where neither gfortran nor
# the C++ compiler contracts to FMA. Elsewhere minqa's Fortran may be
# contracted differently, so only numeric agreement is required.
bitwise <- R.version$arch %in% c("x86_64", "amd64")

# `wellPosed = FALSE` marks runs (e.g. a bound range of 1e-7) whose end
# point is sensitive to rounding, so contracted and uncontracted builds
# legitimately stop at different points.
expectSameAsMinqa <- function(ours, theirs, wellPosed = TRUE) {
  expect_identical(class(ours), class(theirs))
  if (bitwise) {
    expect_identical(ours$par, theirs$par)
    expect_identical(ours$fval, theirs$fval)
    expect_identical(ours$feval, theirs$feval)
    expect_identical(ours$ierr, theirs$ierr)
    expect_identical(ours$msg, theirs$msg)
  } else if (wellPosed) {
    expect_equal(ours$par, theirs$par, tolerance = 1e-4)
  } else {
    expect_identical(ours$ierr, theirs$ierr)
  }
}

test_that("bobyqa() matches minqa::bobyqa() exactly", {
  skip_if_not_installed("minqa")
  cases <- list(
    list(par = c(1, 2), fn = fr, lower = c(0, 0), upper = c(4, 4)),
    list(par = c(-1.2, 1), fn = fr, lower = -Inf, upper = Inf),
    list(par = rep(0.5, 6), fn = extRosen, lower = -2, upper = 2,
         control = list(npt = 13)),
    list(par = seq(0.1, 0.9, length.out = 6), fn = chebyquad, lower = 0,
         upper = 1),
    list(par = c(1, 2), fn = fr, lower = c(0, 0), upper = c(4, 4),
         control = list(maxfun = 50)),
    list(par = c(4, 4), fn = fr, lower = c(0, 3.9999999), upper = c(4, 4),
         wellPosed = FALSE)
  )
  for (cs in cases) {
    ctrl <- if (is.null(cs$control)) list() else cs$control
    ours <- suppressWarnings(bobyqa(cs$par, cs$fn, cs$lower, cs$upper, ctrl))
    theirs <- suppressWarnings(minqa::bobyqa(cs$par, cs$fn, cs$lower,
                                             cs$upper, ctrl))
    expectSameAsMinqa(ours, theirs, wellPosed = !isFALSE(cs$wellPosed))
  }
})

test_that("newuoa() matches minqa::newuoa() exactly", {
  skip_if_not_installed("minqa")
  cases <- list(
    list(par = c(-1.2, 1), fn = fr),
    list(par = rep(0.5, 6), fn = extRosen, control = list(npt = 13)),
    list(par = seq(0.1, 0.9, length.out = 6), fn = chebyquad,
         control = list(npt = 28)),
    list(par = rep(pi, 4), fn = function(x) -(10 - sum(x * seq_along(x))^2),
         control = list(maxfun = 25))
  )
  for (cs in cases) {
    ctrl <- if (is.null(cs$control)) list() else cs$control
    ours <- suppressWarnings(newuoa(cs$par, cs$fn, ctrl))
    theirs <- suppressWarnings(minqa::newuoa(cs$par, cs$fn, ctrl))
    expectSameAsMinqa(ours, theirs)
  }
})

test_that("minqa iprint output matches minqa", {
  skip_if_not_installed("minqa")
  skip_if_not(bitwise)
  # On Windows, Rprintf() (used by minqa) and the C++ snprintf() used by the
  # port's print callback format "%#14.8g" differently (e.g. "196.00" vs
  # "196.00000"), so the text can only be compared elsewhere.
  skip_on_os("windows")
  ours <- capture.output(newuoa(c(1, 2), fr, control = list(iprint = 3)))
  theirs <- capture.output(minqa::newuoa(c(1, 2), fr,
                                         control = list(iprint = 3)))
  expect_identical(ours, theirs)
  ours <- capture.output(bobyqa(c(1, 2), fr, 0, 4, control = list(iprint = 2)))
  theirs <- capture.output(minqa::bobyqa(c(1, 2), fr, 0, 4,
                                         control = list(iprint = 2)))
  expect_identical(ours, theirs)
})

test_that("minqa-style error codes and messages", {
  r <- bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4),
              control = list(maxfun = 50))
  expect_equal(r$ierr, 1)
  expect_match(r$msg, "maximum number of function evaluations")
  expect_s3_class(r, c("bobyqa", "minqa"))
  expect_output(print(r), "parameter estimates")
  r <- newuoa(c(-1.2, 1), fr)
  expect_equal(r$ierr, 0)
  expect_equal(r$msg, "Normal exit from newuoa")
  expect_equal(r$par, c(1, 1), tolerance = 1e-4)
})

test_that("R errors inside fn propagate unchanged", {
  bad <- function(x) if (x[1] > 1.1) stop("boom") else fr(x)
  expect_error(newuoa(c(1, 2), bad), "boom")
  expect_error(bobyqa(c(1, 2), bad, lower = 0, upper = 4), "boom")
  expect_error(bobyqa(c(5, 2), fr, lower = 0, upper = 4),
               "Starting values violate bounds")
})

test_that("extra arguments and names are passed to fn", {
  seen <- NULL
  f <- function(x, a) {
    seen <<- names(x)
    sum((x - a)^2)
  }
  r <- newuoa(c(p = 0.5, q = 0.5), f, a = c(1, 2))
  expect_equal(seen, c("p", "q"))
  expect_equal(r$par, c(1, 2), tolerance = 1e-5)
})

test_that("bobyqa/newuoa C cores give identical results across OpenMP threads", {
  for (solver in c("bobyqa", "newuoa")) {
    res <- RcppTrust:::solver_openmp_stress_test(solver, 6L, 64L, 2L)
    expect_true(res$identical)
  }
})

test_that("NEWUOA's BIGDEN safeguard path matches minqa", {
  skip_if_not_installed("minqa")
  # An ill-conditioned quadratic with npt = (n+1)(n+2)/2 drives NEWUOA into
  # BIGDEN (the denominator-maximizing alternative step).
  dq <- function(x) sum(10^(6 * (seq_along(x) - 1) / (length(x) - 1)) * x^2)
  ctrl <- list(npt = 36, rhobeg = 0.5, rhoend = 1e-10, maxfun = 20000)
  ours <- suppressWarnings(newuoa(rep(1, 7), dq, control = ctrl))
  theirs <- suppressWarnings(minqa::newuoa(rep(1, 7), dq, control = ctrl))
  expectSameAsMinqa(ours, theirs)
})

test_that("BOBYQA's RESCUE path runs deterministically", {
  # A Hilbert-matrix quadratic with a tiny rhoend enters RESCUE, where minqa
  # 1.2.8 itself evaluates uninitialized memory, so there is no minqa
  # reference: the port must converge and repeat exactly.
  hilb <- function(x) {
    n <- length(x)
    drop(x %*% (1 / (outer(0:(n - 1), 0:(n - 1), "+") + 1)) %*% x)
  }
  ctrl <- list(npt = 36, rhobeg = 0.5, rhoend = 1e-10, maxfun = 20000)
  a <- suppressWarnings(bobyqa(rep(1, 7), hilb, -2, 2, control = ctrl))
  b <- suppressWarnings(bobyqa(rep(1, 7), hilb, -2, 2, control = ctrl))
  expect_identical(a, b)
  expect_lt(a$fval, 1e-8)
})

test_that("minqa R wrapper argument handling", {
  expect_warning(newuoa(c(1, 2), fr, control = list(foo = 1)),
                 "unused control arguments")
  expect_warning(newuoa(c(1, 2), fr, control = list(npt = 6)),
                 "not recommended")
  expect_error(newuoa(c(1, 2), function(x) c(1, 2)), "single numeric")
  expect_error(bobyqa(c(1, 2), function(x) "a", 0, 4), "single numeric")
  expect_warning(bobyqa(c(5, 2), fr, lower = 0, upper = 4,
                        control = list(obstop = FALSE)),
                 "adjusted to nearest bound")
  expect_output(bobyqa(c(1, 2), fr, 0, 4, control = list(iprint = 7)),
                "At return")
  r <- newuoa(c(1, 1, 1), function(x) sum(x), control = list(npt = 5))
  expect_equal(r$ierr, 3)
  expect_match(r$msg, "failed to reduce q")
})

test_that("minqa C entry points reject bad input and report callback failures", {
  codes <- RcppTrust:::solver_c_api_edge_test()
  expect_equal(codes[["newuoa_null_opts"]], 0L)
  expect_equal(codes[["bobyqa_null_opts"]], 0L)
  expect_equal(codes[["newuoa_zero_start"]], -3L)
  expect_equal(codes[["bobyqa_negative_n"]], -3L)
  expect_equal(codes[["newuoa_null_result"]], -3L)
  expect_equal(codes[["newuoa_bad_npt"]], 2L)
  expect_equal(codes[["bobyqa_bad_npt"]], 2L)
  expect_equal(codes[["bobyqa_range"]], 4L)
  expect_equal(codes[["newuoa_objfun_error"]], -1L)
  expect_equal(codes[["bobyqa_objfun_error"]], -1L)
  expect_true(codes[["newuoa_nonfinite_x"]] %in% c(-2L, 3L))
})
