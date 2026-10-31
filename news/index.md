# Changelog

## RcppTrust 0.2.1

CRAN release: 2026-10-31

- Errors raised inside the R callbacks of
  [`bobyqa()`](../reference/bobyqa.md),
  [`newuoa()`](../reference/newuoa.md) and
  [`steihaug()`](../reference/steihaug.md) are again reported unchanged
  on clang/libc++ builds (CRAN’s fedora-clang check showed “c++
  exception (unknown reason)”). The error is now stored by value and
  re-thrown afresh instead of going through `std::exception_ptr`
  ([\#3](https://github.com/nlmixr2/RcppTrust/issues/3)).

- The OpenMP thread-safety tests use at most 2 threads and respect
  `OMP_THREAD_LIMIT`, so CRAN no longer sees “OMP: Warning
  [\#96](https://github.com/nlmixr2/RcppTrust/issues/96): Cannot form a
  team” ([\#3](https://github.com/nlmixr2/RcppTrust/issues/3)).

- The
  [`bobyqa()`](../reference/bobyqa.md)/[`newuoa()`](../reference/newuoa.md)
  vs. ‘minqa’ comparison tests compare the permutation-invariant
  chebyquad solutions as a set on platforms without bitwise agreement,
  fixing a test failure on CRAN’s linux-arm64 check where the
  FMA-contracted ‘minqa’ build reached the same minimum with two
  coordinates swapped.

- The NEWUOA `ierr = 3` test is skipped when R lacks long double
  support, since its rounding path depends on
  [`sum()`](https://rdrr.io/r/base/sum.html) accumulating in long double
  (CRAN’s noLD check).

## RcppTrust 0.2.0

CRAN release: 2026-10-27

- The package license is now GPL-2, because the BOBYQA/NEWUOA ports and
  their R wrappers derive from the GPL-2 package ‘minqa’.

- M. J. D. Powell is added as an author for the BOBYQA/NEWUOA algorithms
  and original Fortran code; sources and notices for all ported code
  (‘trust’, ‘basin’, ‘minqa’, Powell’s Fortran) are listed in
  `inst/COPYRIGHTS`.

- New [`steihaug()`](../reference/steihaug.md): a Steihaug truncated
  conjugate-gradient trust-region Newton minimizer ported from the Rust
  crate ‘basin’, with exact-Hessian and matrix-free (Hessian-vector
  product) modes, and a thread-safe C entry point `steihaug_solve_c()`.

- New [`bobyqa()`](../reference/bobyqa.md) and
  [`newuoa()`](../reference/newuoa.md): thread-safe C++ ports of
  Powell’s derivative-free BOBYQA and NEWUOA from ‘minqa’ 1.2.8, drop-in
  replacements for
  [`minqa::bobyqa()`](https://rdrr.io/pkg/minqa/man/bobyqa.html)/[`minqa::newuoa()`](https://rdrr.io/pkg/minqa/man/newuoa.html),
  with C entry points `bobyqa_solve_c()`/`newuoa_solve_c()`. BOBYQA’s
  `RESCUE` step evaluates the intended point rather than minqa’s
  uninitialized one.

- [`.RcppTrustPtr()`](../reference/dot-RcppTrustPtr.md) and
  `RcppTrust.h` gain appended slots for the new C entry points and their
  result-free functions.

- New vignette
  [`vignette("trust-region-methods")`](../articles/trust-region-methods.md).

## RcppTrust 0.1.0

CRAN release: 2026-09-07

- Original release of RcppTRust only containing the port of
  trust::trust()
