# Empirical thread-safety check for the C-objfun path (trust_solve_c()):
# solves many independent Rosenbrock problems once sequentially and once
# with each solve dispatched to a different OpenMP thread, and checks the
# two sets of results are identical. This is what makes trust_solve_c()
# usable from nlmixr2est's parallel C++ loops -- it never touches R's API
# (not thread-safe) and keeps no shared mutable state across calls.

test_that("trust_solve_c() gives identical results run sequentially vs. across OpenMP threads", {
  res <- RcppTrust:::trust_openmp_stress_test(nStarts = 64L, nThreadsRequested = 8L)

  expect_equal(res$maxArgDiff, 0)
  expect_equal(res$maxValDiff, 0)
  expect_true(res$allConvergedMatch)
  expect_true(res$allIterMatch)

  if (res$openmpAvailable) {
    testthat::skip_on_cran()
    if (res$distinctThreadsUsed <= 1) {
      testthat::skip("OpenMP available but only one thread was used by the runtime")
    }
    expect_gt(res$distinctThreadsUsed, 1L)
  } else {
    testthat::skip("Built without OpenMP support -- ran the same check single-threaded")
  }
})
