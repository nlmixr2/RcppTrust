# Ported from upstream trust's tests/hero.R: three scenarios where
# objfun() itself raises an error (as opposed to signaling infeasibility
# via +-Inf), each caught gracefully via the internal .trustTryEval()
# helper and reported through $error the same way upstream's own
# try(objfun(...)) does. Structural values (converged/iterations/
# argument/blather array shapes) are checked against
# tests/hero.Rout.save; $error is only checked for class, since the
# R-shim wraps objfun in an extra closure so the exact call text in the
# condition message differs cosmetically from upstream (e.g. "Error in
# objfun1(theta) : ..." vs "Error in objfun(theta, ...) : ...").

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

test_that("error on the very first call to objfun is caught gracefully", {
  objfun <- function(x) stop("bogus error in objfun")

  expect_warning(
    tout <- trust(objfun, c(3, 1), 1, 5, blather = TRUE),
    "error in first call to objfun"
  )

  expect_true(inherits(tout$error, "try-error"))
  expect_equal(tout$argument, c(3, 1))
  expect_false(tout$converged)
  expect_equal(tout$iterations, 0L)
  expect_null(tout$argpath)  # no blather at all on a first-call error
})

test_that("error mid-loop (on a trial point) is caught gracefully", {
  kiter <- 0
  objfun <- function(x) {
    kiter <<- kiter + 1
    if (kiter == 5) stop("kiter reached 5")
    rosenbrock_objfun(x)
  }

  tout <- trust(objfun, c(3, 1), 1, 5, blather = TRUE)

  expect_true(inherits(tout$error, "try-error"))
  expect_equal(tout$argument, c(1.0323820, 0.1634548), tolerance = 1e-4)
  expect_false(tout$converged)
  expect_equal(tout$iterations, 4L)

  # asymmetric blather lengths: the "pre" arrays recorded one more entry
  # than the "try" arrays, because the failing evaluation happened before
  # its try-array entries could be appended.
  expect_equal(nrow(tout$argpath), 4L)
  expect_equal(nrow(tout$argtry), 3L)
  expect_equal(length(tout$r), 4L)
  expect_equal(length(tout$rho), 3L)
  expect_equal(tout$steptype, c("easy-easy", "easy-easy", "Newton"))
  expect_true(all(tout$accept))
})

test_that("error only on the final post-loop re-evaluation is caught gracefully", {
  kiter <- 0
  objfun <- function(x) {
    kiter <<- kiter + 1
    if (kiter == 5) stop("kiter reached 5")
    rosenbrock_objfun(x)
  }

  expect_warning(
    tout <- trust(objfun, c(3, 1), 1, 5, blather = TRUE, iterlim = 3),
    "error in last call to objfun"
  )

  expect_true(inherits(tout$error, "try-error"))
  expect_equal(tout$argument, c(1.982307, 3.929371), tolerance = 1e-4)
  expect_false(tout$converged)
  expect_equal(tout$iterations, 3L)
  # loop completed normally here, so path/try arrays are the same length
  expect_equal(nrow(tout$argpath), 3L)
  expect_equal(nrow(tout$argtry), 3L)
})
