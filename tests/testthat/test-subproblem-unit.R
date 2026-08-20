# Isolated numeric checks of the ported trust-region subproblem solver
# and its replacement root-finder (boost::math::tools::toms748_solve in
# place of stats::uniroot()), independent of the full trust() loop.
# Exercises internal (non-exported) C++ test shims directly.

test_that("check.objfun.output numeric port matches R semantics", {
  chk <- RcppTrust:::trust_check_numeric_output_test

  ok <- chk(1.5, c(1, 2), c(1, 0, 0, 1), 2, TRUE)
  expect_true(ok$ok)

  r1 <- chk(NaN, NULL, NULL, 2, TRUE)
  expect_false(r1$ok)
  expect_match(r1$message, "NA or NaN")

  r2 <- chk(-Inf, NULL, NULL, 2, TRUE)
  expect_false(r2$ok)
  expect_match(r2$message, "-Inf value in minimization")

  r3 <- chk(Inf, NULL, NULL, 2, FALSE)
  expect_false(r3$ok)
  expect_match(r3$message, "\\+Inf value in maximization")

  # legitimate infeasibility signal: no gradient/hessian required
  r4 <- chk(Inf, NULL, NULL, 2, TRUE)
  expect_true(r4$ok)

  r5 <- chk(1.5, NULL, c(1, 0, 0, 1), 2, TRUE)
  expect_false(r5$ok)
  expect_match(r5$message, "gradient")

  r6 <- chk(1.5, c(1, 2), c(1, 0, 0, NaN), 2, TRUE)
  expect_false(r6$ok)
  expect_match(r6$message, "hessian not having all elements finite")
})

test_that("root-finder (toms748) agrees with stats::uniroot on the same fred()", {
  solve_cpp <- RcppTrust:::trust_solve_lagrange_multiplier_test

  r_solve <- function(beta, gq, r) {
    imin <- beta == 0
    C1 <- sum((gq / beta)[!imin]^2)
    C2 <- sum(gq[imin]^2)
    C3 <- sum(gq^2)
    fred <- function(beep) {
      if (beep == 0) {
        if (C2 > 0) return(-1 / r)
        return(sqrt(1 / C1) - 1 / r)
      }
      sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r
    }
    beta.dn <- sqrt(C2) / r
    beta.up <- sqrt(C3) / r
    if (fred(beta.up) <= 0) {
      beta.up
    } else if (fred(beta.dn) >= 0) {
      beta.dn
    } else {
      stats::uniroot(fred, c(beta.dn, beta.up))$root
    }
  }

  cases <- list(
    list(beta = c(0.5, 2.3, 5.1), gq = c(1.2, -0.7, 0.3), r = 1),
    list(beta = c(0, 1.5, 3.2), gq = c(0.9, -0.4, 0.2), r = 0.8),
    list(beta = c(0.1, 0.2, 0.3), gq = c(1, 1, 1), r = 100)
  )
  set.seed(20260819)
  for (i in 1:5) {
    beta <- sort(runif(5, 0, 10))
    beta[1] <- 0
    cases[[length(cases) + 1]] <- list(beta = beta, gq = rnorm(5), r = runif(1, 0.3, 3))
  }

  for (cs in cases) {
    imin <- cs$beta == 0
    C1 <- sum((cs$gq / cs$beta)[!imin]^2)
    C2 <- sum(cs$gq[imin]^2)
    C3 <- sum(cs$gq^2)
    r_root <- r_solve(cs$beta, cs$gq, cs$r)
    cpp_root <- solve_cpp(cs$beta, cs$gq, cs$r, C1, C2, C3)
    expect_equal(cpp_root, r_root, tolerance = 1e-3)
  }
})

test_that("trust-region subproblem solver matches a direct R port across all step types", {
  solve_cpp <- RcppTrust:::trust_solve_subproblem_test

  r_subproblem <- function(eigval, eigvec, gq, r) {
    is.newton <- FALSE
    if (all(eigval > 0)) {
      ptry <- as.numeric(-eigvec %*% (gq / eigval))
      if (sqrt(sum(ptry^2)) <= r) is.newton <- TRUE
    }
    if (!is.newton) {
      lambda.min <- min(eigval)
      beta <- eigval - lambda.min
      imin <- beta == 0
      C1 <- sum((gq / beta)[!imin]^2)
      C2 <- sum(gq[imin]^2)
      C3 <- sum(gq^2)
      if (C2 > 0 || C1 > r^2) {
        is.hard <- (C2 == 0)
        beta.dn <- sqrt(C2) / r
        beta.up <- sqrt(C3) / r
        fred <- function(beep) {
          if (beep == 0) {
            if (C2 > 0) return(-1 / r)
            return(sqrt(1 / C1) - 1 / r)
          }
          sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r
        }
        root <- if (fred(beta.up) <= 0) {
          beta.up
        } else if (fred(beta.dn) >= 0) {
          beta.dn
        } else {
          stats::uniroot(fred, c(beta.dn, beta.up))$root
        }
        wtry <- gq / (beta + root)
        ptry <- as.numeric(-eigvec %*% wtry)
        steptype <- if (is.hard) "hard-easy" else "easy-easy"
      } else {
        wtry <- gq / beta
        wtry[imin] <- 0
        ptry <- as.numeric(-eigvec %*% wtry)
        utry <- sqrt(r^2 - sum(ptry^2))
        if (utry > 0) {
          vtry <- eigvec[, imin, drop = FALSE][, 1]
          ptry <- ptry + utry * vtry
        }
        steptype <- "hard-hard"
      }
    } else {
      steptype <- "Newton"
    }
    list(ptry = ptry, steptype = steptype)
  }

  cases <- list(
    list(B = diag(c(4, 2)), g = c(1, 1), r = 5, label = "Newton"),
    list(B = diag(c(2, -2)), g = c(0.5, 0.5), r = 1, label = "easy-easy"),
    list(B = diag(c(2, -2)), g = c(0, 0), r = 1, label = "hard-hard")
  )
  for (cs in cases) {
    eout <- eigen(cs$B, symmetric = TRUE)
    gq <- as.numeric(t(eout$vectors) %*% cs$g)
    rres <- r_subproblem(eout$values, eout$vectors, gq, cs$r)
    cres <- solve_cpp(eout$values, eout$vectors, gq, cs$r)
    expect_equal(cres$steptype, rres$steptype, info = cs$label)
    expect_equal(cres$ptry, rres$ptry, tolerance = 1e-3, info = cs$label)
  }
})
