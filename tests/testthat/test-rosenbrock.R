# Ported from upstream trust's tests/foo.R; golden values below are
# transcribed from tests/foo.Rout.save. Root-find (uniroot -> toms748)
# differences accumulate over iterations, so numeric comparisons use a
# loose tolerance; structural fields (steptype, accept, iterations,
# converged) are checked exactly.

rosenbrock_objfun <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 2)
  f <- expression(100 * (x2 - x1^2)^2 + (1 - x1)^2)
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

test_that("Rosenbrock example matches upstream trust::trust()", {
  tout <- trust(rosenbrock_objfun, c(3, 1), 1, 5, blather = TRUE)

  expect_true(tout$converged)
  expect_equal(tout$iterations, 21L)
  expect_equal(tout$argument, c(1, 1), tolerance = 1e-4)
  expect_equal(tout$value, 5.160801e-15, tolerance = 1e-4)

  expect_equal(
    tout$steptype,
    c(
      "easy-easy", "easy-easy", "Newton", "Newton", "easy-easy", "easy-easy",
      "easy-easy", "Newton", "Newton", "Newton", "easy-easy", "easy-easy",
      "Newton", "Newton", "easy-easy", "easy-easy", "Newton", "Newton",
      "Newton", "Newton", "Newton"
    )
  )
  expect_equal(
    tout$accept,
    c(
      TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
      TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    )
  )
  expect_equal(tout$r, c(
    1, 2, 4, 4, 1, 0.25, 0.5, 1, 1, 1, 0.25, 0.5, 0.5, 0.5, 0.125, 0.25, 0.5,
    0.5, 0.5, 0.5, 0.5
  ))

  # Constrained non-Newton accepted steps sit exactly on the trust region
  # boundary.
  ratio <- (tout$stepnorm / tout$r)[tout$accept & tout$steptype != "Newton"]
  expect_equal(ratio, rep(1, length(ratio)), tolerance = 1e-3)
})
