# Ported from upstream trust's tests/bar.R and tests/baz.R. Both start at
# a point where the trust region subproblem hits the "hard-hard" case, in
# which the algorithm adds a multiple of the eigenvector spanning the
# minimal-eigenvalue eigenspace to reach the trust region boundary. That
# eigenvector is only defined up to sign (an ordinary eigendecomposition
# ambiguity at exact eigenvalue ties, independent of the algorithm), and
# RcppTrust uses RcppArmadillo/LAPACK's eig_sym rather than R's eigen(),
# which can pick the opposite sign than upstream did -- confirmed by
# direct comparison against trust::trust() to be a pure sign flip of the
# degenerate component (identical objective value, steptype sequence,
# iteration count). So argument comparisons here go through
# expect_equal_up_to_sign() on the component known to be affected.

bar_objfun <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 2)
  f <- expression(x1^2 - x2^2 + (1 / 100) * (x1^4 + x2^4))
  g1 <- D(f, "x1")
  g2 <- D(f, "x2")
  h11 <- D(g1, "x1")
  h12 <- D(g1, "x2")
  h22 <- D(g2, "x2")
  x1 <- x[1]
  x2 <- x[2]
  f <- eval(f)
  g <- c(eval(g1), eval(g2))
  B <- rbind(c(eval(h11), eval(h12)), c(eval(h12), eval(h22)))
  list(value = f, gradient = g, hessian = B)
}

baz_objfun <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 2)
  f <- expression((1 / 10) * x1 + x1^2 - x2^2 + (1 / 100) * (x1^4 + x2^4))
  g1 <- D(f, "x1")
  g2 <- D(f, "x2")
  h11 <- D(g1, "x1")
  h12 <- D(g1, "x2")
  h22 <- D(g2, "x2")
  x1 <- x[1]
  x2 <- x[2]
  f <- eval(f)
  g <- c(eval(g1), eval(g2))
  B <- rbind(c(eval(h11), eval(h12)), c(eval(h12), eval(h22)))
  list(value = f, gradient = g, hessian = B)
}

test_that("bar: hard-hard case at a saddle point", {
  tout <- trust(bar_objfun, c(0, 0), 1, 5, blather = TRUE)

  expect_true(tout$converged)
  expect_equal(tout$iterations, 6L)
  expect_equal(tout$value, -25, tolerance = 1e-6)
  expect_equal(tout$argument[1], 0, tolerance = 1e-4)
  expect_equal(abs(tout$argument[2]), 7.071068, tolerance = 1e-4)
  expect_equal(
    tout$steptype,
    c("hard-hard", "easy-easy", "easy-easy", "Newton", "Newton", "Newton")
  )
  expect_true(all(tout$accept))
  expect_equal(tout$r, c(1, 2, 4, 4, 4, 4))
})

test_that("baz: hard-hard case broken slightly asymmetric by a linear term", {
  tout <- trust(baz_objfun, c(0, 0), 1, 5, blather = TRUE)

  expect_true(tout$converged)
  expect_equal(tout$iterations, 6L)
  expect_equal(tout$value, -25.0025, tolerance = 1e-5)
  expect_equal(tout$argument[1], -0.0499975, tolerance = 1e-4)
  expect_equal(abs(tout$argument[2]), 7.0710678, tolerance = 1e-4)
  expect_equal(
    tout$steptype,
    c("hard-hard", "easy-easy", "easy-easy", "Newton", "Newton", "Newton")
  )
  expect_true(all(tout$accept))
})
