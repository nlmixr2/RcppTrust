# Ported from upstream trust's tests/poo.R: verifies that the built-in
# `parscale` argument is equivalent to manually reparameterizing the
# objective function -- trust(objfun, parinit, ..., parscale = parscale)
# should behave identically to calling trust() on
# function(x) objfun(x / parscale - shift) starting from
# parscale * (parinit + shift).

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

test_that("parscale argument matches manual reparameterization", {
  parinit <- c(3, 1)

  out <- trust(rosenbrock_objfun, parinit, 1, 1e5, blather = TRUE)
  expect_true(out$converged)
  expect_equal(length(out$r), 21L)

  parscale <- c(5, 1)
  shift <- 4
  theta <- parscale * (parinit + shift)

  pobjfun <- function(x) {
    o <- rosenbrock_objfun(x / parscale - shift)
    o$gradient <- o$gradient / parscale
    o$hessian <- o$hessian / outer(parscale, parscale)
    o
  }

  pout <- trust(pobjfun, theta, 1, 1e5, blather = TRUE)
  expect_true(pout$converged)
  expect_equal(length(pout$r), 25L)
  expect_equal(out$argument, pout$argument / parscale - shift, tolerance = 1e-4)

  qout <- trust(rosenbrock_objfun, parinit, 1, 1e5, parscale = parscale, blather = TRUE)
  expect_true(qout$converged)
  expect_equal(length(qout$r), 25L)

  expect_equal(pout$valpath, qout$valpath, tolerance = 1e-4)
  transpath <- sweep(sweep(pout$argpath, 2, parscale, "/"), 2, shift)
  expect_equal(unname(transpath), unname(qout$argpath), tolerance = 1e-4)
})
