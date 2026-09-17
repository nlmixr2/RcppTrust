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

test_that("every Steihaug stopping rule is reachable from R", {
  run <- function(...) steihaug(c(-1.2, 1), fr, grr, hr,
                                control = c(list(gradTol = NULL), list(...)))
  expect_equal(run(relGradTol = 1e-6)$message, "RelativeGradientTolerance")
  expect_equal(run(stepTol = 1e-8)$message, "ParamTolerance")
  expect_equal(run(relStepTol = 1e-8)$message, "RelativeParamTolerance")
  expect_equal(run(costTol = 1e-14)$message, "CostTolerance")
  expect_equal(run(relCostTol = 0.5)$message, "RelativeCostTolerance")
  expect_equal(run(maxCostEvals = 10)$message, "MaxCostEvals")
  expect_equal(run(maxGradEvals = 5)$message, "MaxGradientEvals")
  r <- run(maxit = 1000)
  expect_equal(r$message, "SolverConverged")
  expect_true(r$converged)
  # General forcing exponent (not 0, 0.5 or 1) and a CG cap.
  r <- steihaug(c(-1.2, 1), fr, grr, hr,
                control = list(theta = 0.7, kappa = 0.1, cgMaxit = 1,
                               maxit = 50))
  expect_equal(r$message, "MaxIter")
  r <- steihaug(c(-1.2, 1), fr, grr, hr, control = list(theta = 0))
  expect_true(r$converged)
  # A NaN gradient is a failure, not convergence.
  r <- steihaug(c(1, 1), fr, function(x) c(NaN, 0), hr)
  expect_equal(r$message, "SolverFailed")
  expect_false(r$converged)
})

test_that("steihaug() validates its inputs", {
  expect_error(steihaug(c(1, NA), fr, grr, hr), "finite")
  expect_error(steihaug(c(1, 1), 1, grr, hr), "fn must")
  expect_error(steihaug(c(1, 1), fr, 1, hr), "gr must")
  expect_error(steihaug(c(1, 1), fr, grr, 1), "hess must")
  expect_error(steihaug(c(1, 1), fr, grr, hessvec = 1), "hessvec must")
  expect_error(steihaug(c(1, 1), fr, grr, hr, control = 1), "list")
  expect_error(steihaug(c(1, 1), fr, grr, hr, control = list(rinit = NA)),
               "rinit")
  expect_error(steihaug(c(1, 1), fr, grr, hr, control = list(gradTol = "a")),
               "gradTol")
  expect_error(steihaug(c(1, 1), function(x) c(1, 2), grr, hr), "single")
  expect_error(steihaug(c(-1.2, 1), fr, grr, hessvec = function(x, v) 1),
               "length 2")
  expect_error(steihaug(c(-1.2, 1), fr, grr,
                        hessvec = function(x, v) stop("hvx")),
               "hvx")
})

test_that("Steihaug C entry point rejects bad input and reports failures", {
  codes <- RcppTrust:::solver_c_api_edge_test()
  expect_equal(codes[["steihaug_ok"]], 4L)
  expect_equal(codes[["steihaug_null_opts"]], -3L)
  expect_equal(codes[["steihaug_bad_eta"]], -3L)
  expect_equal(codes[["steihaug_nan_tol"]], -3L)
  expect_equal(codes[["steihaug_objfun_error"]], -1L)
  expect_equal(codes[["steihaug_hessvec_error"]], -2L)
  expect_equal(codes[["steihaug_zero_n"]], -3L)
})
