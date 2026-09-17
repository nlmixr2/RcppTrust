# RcppTrust 0.2.0

* The package license is now GPL-2, because the new BOBYQA/NEWUOA ports
  derive from the GPL-2 package 'minqa'. Notices for the 'trust' and
  'basin' code are retained in `inst/COPYRIGHTS`.

* New `steihaug()`: a Steihaug truncated conjugate-gradient trust-region
  Newton minimizer ported from the Rust crate 'basin', with exact-Hessian
  and matrix-free (Hessian-vector product) modes, and a thread-safe C
  entry point `steihaug_solve_c()`.

* New `bobyqa()` and `newuoa()`: thread-safe C++ ports of Powell's
  derivative-free BOBYQA and NEWUOA from 'minqa' 1.2.8, drop-in
  replacements for `minqa::bobyqa()`/`minqa::newuoa()`, with C entry points
  `bobyqa_solve_c()`/`newuoa_solve_c()`. BOBYQA's `RESCUE` step evaluates
  the intended point rather than minqa's uninitialized one.

* `.RcppTrustPtr()` and `RcppTrust.h` gain appended slots for the new C
  entry points and their result-free functions.

* New vignette `vignette("trust-region-methods")`.
