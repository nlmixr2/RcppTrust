## Resubmission (0.2.1)

This is a quick resubmission (3 days after 0.2.0) to fix the check
problems CRAN asked me to correct on
https://cran.r-project.org/web/checks/check_results_RcppTrust.html:

* r-devel-linux-x86_64-fedora-clang ERROR ("c++ exception (unknown
  reason)" instead of the user's R error from `newuoa()`/`steihaug()`
  callbacks). Errors raised in R callbacks are now stored by value and
  re-thrown afresh instead of via `std::exception_ptr`, which did not
  survive the libc++ build.

* "OMP: Warning #96: Cannot form a team with 8 threads": the OpenMP
  thread-safety tests now use at most 2 threads and respect
  `OMP_THREAD_LIMIT`.

* linux-arm64 additional issue: the `bobyqa()` vs. 'minqa' chebyquad
  comparison now compares the permutation-invariant solution as a set
  where bitwise agreement is not expected (FMA-contracted 'minqa' reached
  the same minimum with two coordinates swapped).

* noLD additional issue: the NEWUOA `ierr = 3` test is skipped when R
  lacks long double support, since that rounding path depends on `sum()`
  accumulating in long double.

These fixes were verified on:

* R-hub `nold` (R-devel, --disable-long-double):
  https://github.com/nlmixr2/RcppTrust/actions/runs/40872715805
  -- Status: OK, [ FAIL 0 | PASS 293 ]
* R-hub `ubuntu-clang` (R-devel, clang 23):
  https://github.com/nlmixr2/RcppTrust/actions/runs/40872932771
  -- Status: OK, [ FAIL 0 | PASS 295 ]
* linux-arm64 using CRAN's `rcheckserver` image (mirrors the CRAN
  linux-arm64 check): tests pass.
* GitHub Actions: macOS (Apple clang / libc++), Windows, and Ubuntu
  R-devel/release/oldrel-1.

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 3 -- this release fixes the CRAN check
  problems listed above.
