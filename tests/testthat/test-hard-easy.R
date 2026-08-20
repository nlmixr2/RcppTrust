# Ported from upstream trust's tests/qux.R. Exercises the "hard-easy"
# case on the first iteration and "hard-hard" on the second. Golden
# values transcribed from tests/qux.Rout.save. See test-hard-hard.R for
# why the degenerate-direction component of argument is compared up to
# sign rather than exactly.

qux_objfun <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 2)
  f <- expression(5 * x1 + x1^2 - x2^2 + (1 / 100) * (x1^4 + x2^4))
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

test_that("qux: hard-easy then hard-hard trust-region subproblem cases", {
  tout <- trust(qux_objfun, c(0, 0), 1, 5, blather = TRUE)

  expect_true(tout$converged)
  expect_equal(tout$iterations, 8L)
  expect_equal(tout$value, -30.93159, tolerance = 1e-4)
  expect_equal(tout$argument[1], -2.266988, tolerance = 1e-3)
  expect_equal(abs(tout$argument[2]), 7.071068, tolerance = 1e-3)

  expect_equal(
    tout$steptype,
    c(
      "hard-easy", "hard-hard", "easy-easy", "Newton", "Newton", "Newton",
      "Newton", "Newton"
    )
  )
  expect_true(all(tout$accept))
  expect_equal(tout$r, c(1, 2, 4, 4, 4, 4, 4, 4))

  ratio <- (tout$stepnorm / tout$r)[tout$accept & tout$steptype != "Newton"]
  expect_equal(ratio, rep(1, length(ratio)), tolerance = 1e-3)
})
