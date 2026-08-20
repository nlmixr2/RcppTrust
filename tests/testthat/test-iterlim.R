# Ported from upstream trust's tests/goo.R: an unbounded objective (no
# minimum exists), verifying iterlim is respected and the algorithm
# degrades gracefully (every step type stays "easy-easy", every step
# lands exactly on the trust region boundary, radius keeps growing to
# rmax) rather than erroring. Golden values from tests/goo.Rout.save.

goo_objfun <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 2)
  f <- expression(x1^2 - x2^2)
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

test_that("goo: iterlim is respected on an unbounded objective", {
  goo <- trust(goo_objfun, c(3, 1), 1, 5, blather = TRUE, iterlim = 20)

  expect_false(goo$converged)
  expect_equal(goo$iterations, 20L)
  expect_equal(length(goo$r), 20L)
  expect_true(all(goo$steptype == "easy-easy"))
  expect_true(all(goo$accept))
  expect_true(all(goo$rho > 0.99))

  # radius grows geometrically to rmax = 5 and stays there
  expect_equal(goo$r[1:3], c(1, 2, 4))
  expect_true(all(goo$r[4:20] == 5))

  ratio <- (goo$stepnorm / goo$r)[goo$accept & goo$steptype != "Newton"]
  expect_equal(ratio, rep(1, length(ratio)), tolerance = 1e-3)
})
