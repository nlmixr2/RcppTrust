fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
grr <- function(x) c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
                     200 * (x[2] - x[1]^2))
hr <- function(x) matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
                           -400 * x[1], 200), 2, 2)
hv <- function(x, v) drop(hr(x) %*% v)

test_that("steihaug() minimizes Rosenbrock with basin's evaluation pattern", {
  r <- steihaug(c(-1.2, 1), fr, grr, hr)
  expect_true(r$converged)
  expect_equal(r$message, "GradientTolerance")
  expect_equal(r$par, c(1, 1), tolerance = 1e-8)
  # Reference counts from basin 1.12.0 (TrustRegion<Steihaug>,
  # absolute gradient tolerance 1e-8): 27 iterations, 31 cost,
  # 28 gradient and 27 Hessian evaluations.
  expect_equal(r$iterations, 27L)
  expect_equal(unname(r$counts), c(31, 28, 27, 0))
})

test_that("matrix-free mode follows the same path as exact mode", {
  e <- steihaug(c(-1.2, 1), fr, grr, hr)
  m <- steihaug(c(-1.2, 1), fr, grr, hessvec = hv)
  expect_identical(e$par, m$par)
  expect_identical(e$iterations, m$iterations)
  expect_equal(unname(m$counts), c(31, 28, 0, 76))
})

test_that("control settings and termination reasons", {
  r <- steihaug(c(-1.2, 1), fr, grr, hr, control = list(maxit = 5))
  expect_equal(r$message, "MaxIter")
  expect_false(r$converged)
  expect_equal(r$iterations, 5L)

  r <- steihaug(c(-1.2, 1), fr, grr, hr,
                control = list(gradTol = NULL, targetCost = 1e-6))
  expect_equal(r$message, "TargetCost")
  expect_lte(r$value, 1e-6)

  r <- steihaug(c(-1.2, 1), fr, grr, hr, control = list(trace = TRUE))
  expect_s3_class(r$trace, "data.frame")
  expect_equal(sum(r$trace$accept), r$iterations)

  expect_error(steihaug(c(-1.2, 1), fr, grr, hr, control = list(eta = 0.5)),
               "invalid")
  expect_warning(steihaug(c(-1.2, 1), fr, grr, hr, control = list(foo = 1)),
                 "unused control")
  expect_error(steihaug(c(-1.2, 1), fr, grr), "hess")
  expect_error(steihaug(c(-1.2, 1), fr, grr, hr, control = list(maxit = Inf)),
               "maxit")
  expect_error(steihaug(numeric(0), fr, grr, hr), "non-empty")
})

test_that("errors and malformed returns in R callbacks propagate", {
  expect_error(steihaug(c(-1.2, 1), fr, function(x) stop("gradfail"), hr),
               "gradfail")
  expect_error(steihaug(c(-1.2, 1), fr, function(x) 1, hr), "length 2")
  expect_error(steihaug(c(-1.2, 1), fr, grr, function(x) diag(3)), "2 by 2")
})

test_that("infinite trial values are rejected, not accepted", {
  f <- function(x) if (any(x <= 0)) Inf else sum(x - log(x))
  g <- function(x) 1 - 1 / x
  h <- function(x) diag(1 / x^2, length(x))
  r <- steihaug(c(3, 0.2), f, g, h)
  expect_true(r$converged)
  expect_equal(r$par, c(1, 1), tolerance = 1e-6)
})

test_that("steihaug C core gives identical results across OpenMP threads", {
  for (solver in c("steihaug", "steihaug_hessvec")) {
    res <- RcppTrust:::solver_openmp_stress_test(solver, 8L, 64L, 2L)
    expect_true(res$identical)
    expect_gt(res$nNearOptimum, 0L)
  }
})
