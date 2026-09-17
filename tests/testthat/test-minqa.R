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

expectSameAsMinqa <- function(ours, theirs) {
  expect_identical(class(ours), class(theirs))
  if (bitwise) {
    expect_identical(ours$par, theirs$par)
    expect_identical(ours$fval, theirs$fval)
    expect_identical(ours$feval, theirs$feval)
    expect_identical(ours$ierr, theirs$ierr)
    expect_identical(ours$msg, theirs$msg)
  } else {
    expect_equal(ours$par, theirs$par, tolerance = 1e-4)
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
    list(par = c(4, 4), fn = fr, lower = c(0, 3.9999999), upper = c(4, 4))
  )
  for (cs in cases) {
    ctrl <- if (is.null(cs$control)) list() else cs$control
    ours <- suppressWarnings(bobyqa(cs$par, cs$fn, cs$lower, cs$upper, ctrl))
    theirs <- suppressWarnings(minqa::bobyqa(cs$par, cs$fn, cs$lower,
                                             cs$upper, ctrl))
    expectSameAsMinqa(ours, theirs)
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
