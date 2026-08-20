# Ported from upstream trust's tests/hoo.R (also the second \examples
# block in trust.Rd): a log-barrier-type objective with domain restricted
# to the open unit ball, exercising the objfun-returns-Inf-for-infeasible
# convention, and the "parinit not feasible" hard-stop for an infeasible
# starting point. Golden values from tests/hoo.Rout.save.

restricted_objfun_factory <- function(mu, d) {
  function(x) {
    stopifnot(is.numeric(x))
    stopifnot(length(x) == d)
    normxsq <- sum(x^2)
    omnormxsq <- 1 - normxsq
    if (normxsq >= 1) return(list(value = Inf))
    f <- sum(x * mu) - log(omnormxsq)
    g <- mu + 2 * x / omnormxsq
    B <- 4 * outer(x, x) / omnormxsq^2 + 2 * diag(d) / omnormxsq
    list(value = f, gradient = g, hessian = B)
  }
}

test_that("hoo: restricted-domain objective converges to the documented solution", {
  d <- 5
  mu <- seq(1, d)
  objfun <- restricted_objfun_factory(mu, d)

  whoop <- trust(objfun, rep(0, d), 1, 100, blather = TRUE)

  expect_true(whoop$converged)
  expect_lt(max(abs(whoop$gradient)), 1e-8)
  expect_equal(length(whoop$r), 11L)
  expect_equal(
    whoop$argument, c(-0.118, -0.236, -0.354, -0.472, -0.589),
    tolerance = 1e-3
  )
  # golden value was printed at options(digits = 3), so it's only accurate
  # to 3 significant figures
  expect_equal(1 - sqrt(sum(whoop$argument^2)), 0.126, tolerance = 5e-3)

  ratio <- (whoop$stepnorm / whoop$r)[whoop$accept & whoop$steptype != "Newton"]
  expect_equal(ratio, rep(1, length(ratio)), tolerance = 1e-3)
})

test_that("hoo: a harder, more ill-conditioned instance still converges", {
  d <- 5
  mu <- 10 * seq(1, d)
  objfun <- restricted_objfun_factory(mu, d)

  whoop <- trust(objfun, rep(0, d), 1, 100, blather = TRUE)

  expect_true(whoop$converged)
  expect_lt(max(abs(whoop$gradient)), 1e-9)
  expect_equal(length(whoop$r), 15L)
  expect_equal(
    whoop$argument, c(-0.133, -0.266, -0.399, -0.532, -0.665),
    tolerance = 3e-3
  )
  # golden value was printed at options(digits = 3), so it's only accurate
  # to 3 significant figures
  expect_equal(1 - sqrt(sum(whoop$argument^2)), 0.0134, tolerance = 1e-2)
})

test_that("hoo: an infeasible starting point is a hard, uncaught error", {
  d <- 5
  mu <- seq(1, d)
  objfun <- restricted_objfun_factory(mu, d)

  expect_error(
    trust(objfun, rep(0.5, d), 1, 100, blather = TRUE),
    "parinit not feasible"
  )
})
