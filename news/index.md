# Changelog

## RcppTrust 0.2.0

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
