# Ported from upstream trust's tests/fred.R and tests/gred.R: fits a
# four-way contingency table Poisson log-linear model (86 parameters) via
# trust() and cross-checks against glm.fit(). Reuses the frozen fred.txt
# data (copied byte-for-byte from upstream, since it was generated with
# sample(), whose algorithm changed across R versions -- unlike rnorm()'s
# inversion method, which has been stable, so theta.true is regenerated
# from set.seed(42) exactly as upstream does).

test_that("fred/gred: Poisson log-linear model fit matches glm.fit()", {
  skip_if_not_installed("stats")

  data <- read.table(testthat::test_path("fred.txt"), header = TRUE)
  x <- data$x
  data$x <- NULL
  m <- as.matrix(data)
  dimnames(m) <- NULL

  set.seed(42)
  theta.true <- 0.25 * rnorm(ncol(m))
  theta.true <- round(theta.true, 5)

  objfun <- function(theta) {
    eta <- as.numeric(m %*% theta)
    p <- exp(eta)
    f <- sum(x * eta - p)
    g <- as.numeric(t(x - p) %*% m)
    B <- sweep(-m, 1, p, "*")
    B <- t(m) %*% B
    list(value = f, gradient = g, hessian = B)
  }

  # verify the analytic gradient/hessian by finite differences, as
  # upstream's own test does, before trusting them to trust()
  sally <- objfun(theta.true)
  epsilon <- 1e-8
  mygrad <- double(length(theta.true))
  for (i in seq_along(mygrad)) {
    theta.eps <- theta.true
    theta.eps[i] <- theta.true[i] + epsilon
    mygrad[i] <- (objfun(theta.eps)$value - sally$value) / epsilon
  }
  expect_equal(sally$gradient, mygrad, tolerance = length(mygrad) * epsilon)

  # fred.R: maximize (minimize = FALSE) matches glm.fit()
  fred <- trust(objfun, theta.true, 1, sqrt(ncol(m)), minimize = FALSE)
  fran <- stats::glm.fit(m, x, family = stats::poisson(), intercept = FALSE)
  expect_equal(unname(fran$coefficients), fred$argument, tolerance = 1e-4)

  fred0 <- trust(objfun, rep(0, length(theta.true)), 1, sqrt(ncol(m)), minimize = FALSE)
  expect_true(fred0$converged)

  fredNeg5 <- trust(objfun, rep(-5, length(theta.true)), 1, sqrt(ncol(m)), minimize = FALSE)
  expect_true(fredNeg5$converged)
})

test_that("gred: forgetting minimize = FALSE fails to converge (unbounded below)", {
  data <- read.table(testthat::test_path("fred.txt"), header = TRUE)
  x <- data$x
  data$x <- NULL
  m <- as.matrix(data)
  dimnames(m) <- NULL

  set.seed(42)
  theta.true <- 0.25 * rnorm(ncol(m))
  theta.true <- round(theta.true, 5)

  objfun <- function(theta) {
    eta <- as.numeric(m %*% theta)
    p <- exp(eta)
    f <- sum(x * eta - p)
    g <- as.numeric(t(x - p) %*% m)
    B <- sweep(-m, 1, p, "*")
    B <- t(m) %*% B
    list(value = f, gradient = g, hessian = B)
  }

  gred <- trust(objfun, theta.true, 1, sqrt(ncol(m)), blather = TRUE)

  expect_false(gred$converged)
  expect_equal(gred$iterations, 100L)
  expect_true(all(gred$steptype == "easy-easy"))
  expect_equal(gred$r[length(gred$r)], 9.27, tolerance = 1e-2)

  ratio <- (gred$stepnorm / gred$r)[gred$accept & gred$steptype != "Newton"]
  expect_equal(ratio, rep(1, length(ratio)), tolerance = 1e-3)
})
