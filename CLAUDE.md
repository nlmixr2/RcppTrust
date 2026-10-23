# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when
working with code in this repository.

## Overview

**RcppTrust** is an R package of thread-safe C++ optimizers built for
use from parallel C++ code in the nlmixr2 ecosystem (e.g. nlmixr2est).
Each optimizer is a C-callable core with no R API calls and no global
state, wrapped by a thin R function:

- `trust()` -- port of `trust::trust()` (Geyer); same arguments, defaults
  and return value.
- `steihaug()` -- Steihaug truncated conjugate-gradient trust-region
  Newton method, ported from the Rust crate 'basin' (exact Hessian or
  matrix-free Hessian-vector product modes).
- `bobyqa()` / `newuoa()` -- line-by-line ports of Powell's
  derivative-free BOBYQA/NEWUOA as shipped in 'minqa' 1.2.8; drop-in
  replacements for `minqa::bobyqa()` / `minqa::newuoa()` (`print.minqa`
  S3 method included).
- `.RcppTrustPtr()` -- positionally-indexed function-pointer table for
  downstream packages (see "Cross-package C API" below).

License is GPL-2 (because of the minqa-derived code). Sources and
notices for all ported code are in `inst/COPYRIGHTS`. Depends on
R >= 4.2.0, C++17; Imports Rcpp; LinkingTo Rcpp, RcppArmadillo, BH.

## Build and Development Commands

```sh
# Install
R CMD INSTALL .

# Regenerate NAMESPACE/man (roxygen2) and RcppExports after editing
# roxygen comments or // [[Rcpp::export]] functions
Rscript -e "Rcpp::compileAttributes(); devtools::document()"

# Run all tests (remove stale objects first)
find src -name "*.so" -o -name "*.o" | xargs rm -f 2>/dev/null; NOT_CRAN=true Rscript -e "devtools::test()"

# Run one test file (name without "test-" prefix and ".R" suffix)
NOT_CRAN=true Rscript -e "devtools::test(filter='minqa')"
NOT_CRAN=true Rscript -e "devtools::test(filter='steihaug')"

# R CMD check
Rscript -e "devtools::check()"

# Rebuild README.md after editing README.Rmd
Rscript -e "devtools::build_readme()"
```

CI (`.github/workflows/`): `R-CMD-check.yaml` (macOS, Windows, Ubuntu
release/devel/oldrel-1), `test-coverage.yaml` (codecov),
`pkgdown.yaml`, `rhub.yaml`.

## Architecture

### Layout per solver

Each solver is split into a thread-safe core, a public C header with
its option/result structs and C entry point, and an Rcpp R interface:

| Solver   | Core (src/)                                   | Public header (inst/include/) | R interface                                   |
|----------|-----------------------------------------------|-------------------------------|-----------------------------------------------|
| trust    | `trust_core.cpp`, `trust_subproblem.cpp`, `trust_rootfind.cpp`, `trust_check.cpp`, `trust_result.cpp` | `trust_types.h` | `trust_r_interface.cpp`, `trust_robjfun.cpp`, `R/trust.R` |
| Steihaug | `steihaug_core.cpp`                           | `steihaug_types.h`            | `steihaug_r_interface.cpp`, `R/steihaug.R`    |
| BOBYQA   | `minqa_bobyqa.cpp`, `minqa_common.cpp`        | `minqa_types.h`               | `minqa_r_interface.cpp`, `R/minqa.R`          |
| NEWUOA   | `minqa_newuoa.cpp`, `minqa_common.cpp`        | `minqa_types.h`               | `minqa_r_interface.cpp`, `R/minqa.R`          |

- `src/*_types.h` are forwarding headers that `#include
  "../inst/include/*_types.h"`; the canonical struct definitions live
  only in `inst/include/` so the package and consumers share one copy.
  Edit the `inst/include/` version.
- C entry points: `trust_solve_c` / `trust_result_free`,
  `steihaug_solve_c` / `steihaug_result_free`, `bobyqa_solve_c` /
  `newuoa_solve_c` / `minqa_result_free`.
- `R/trust.R` keeps `.trustTryEval()` in R deliberately: an error in the
  user's objfun must produce a genuine R "try-error" exactly as upstream
  `trust::trust()` does. Do not reimplement it in C++.
- `R/c-api.R` holds only roxygen docs for the C API (`minqa_c_api` Rd).

### Cross-package C API (do not break the ABI)

- `inst/include/RcppTrust.h` is a header-only, positionally-indexed
  function-pointer registration (same pattern as n1qn1c/lbfgsb3c/rxode2).
  Consumers never link against RcppTrust's shared library; they fill
  pointers at load time from `.RcppTrustPtr()`.
- `src/trust_register.cpp`: `trust_ptr()` builds that list;
  `trust_register_ccallable()` (`// [[Rcpp::init]]`) also registers each
  pointer with `R_RegisterCCallable`. There is no hand-written
  `src/init.c` -- Rcpp generates `R_init_RcppTrust` in `RcppExports.cpp`.
- **Slot order is ABI. Only ever append new slots; never reorder or
  remove.** When adding a C entry point, update in lockstep:
  `trust_ptr()` names/fns arrays, `trust_register_ccallable()`, the
  typedefs/externs/`iniRcppTrustPtrs0()`/`iniRcppTrust` macro in
  `RcppTrust.h`, the slot list in the `R/c-api.R` docs, and
  NEWS.md. Every solver must also expose its matching `*_result_free`.
- Slots beyond an older RcppTrust's list length are left NULL;
  consumers must check pointers before use.

### Thread safety and numerical reproducibility

- Core code must not call the R API (no `Rprintf`, `Rf_error`, R
  allocation) and must not use global/static mutable state.
- The BOBYQA/NEWUOA ports keep Powell's variable names, workspace
  partitioning, 1-based index macros, goto labels and floating-point
  operation order so results are bitwise identical to the Fortran.
  `minqa_bobyqa.cpp`, `minqa_newuoa.cpp` and `steihaug_core.cpp` forbid
  FMA contraction per translation unit; keep that in place and do not
  "tidy" the operation order.
- OpenMP is enabled via `$(SHLIB_OPENMP_CXXFLAGS)` in `src/Makevars`
  and `src/Makevars.win` (keep both in sync).

### Test-only native code

`src/trust_test_objfuns.cpp` and `src/solver_test_objfuns.cpp` provide
native objectives, OpenMP stress tests (`trust_openmp_stress_test`,
`solver_openmp_stress_test`), C API edge tests
(`solver_c_api_edge_test`) and subproblem unit hooks. They are
non-exported and reached from tests via `RcppTrust:::`.

## Generated Files (do not edit manually)

- `src/RcppExports.cpp`, `R/RcppExports.R` -- `Rcpp::compileAttributes()`
- `NAMESPACE`, `man/*.Rd` -- roxygen2
- `README.md` -- from `README.Rmd`

## Testing Conventions

- testthat edition 3; files in `tests/testthat/` named by solver or
  feature (`test-minqa.R`, `test-steihaug.R`, `test-rosenbrock.R`,
  `test-hard-hard.R`, `test-subproblem-unit.R`, `test-c-objfun-threadsafety.R`, ...).
- Behavior is checked against the reference packages: guard comparisons
  with `skip_if_not_installed("trust")` / `skip_if_not_installed("minqa")`
  since these are only in Suggests.
- Thread-safety tests compare sequential vs. OpenMP results for exact
  equality; the multi-thread assertion is behind `skip_on_cran()` and
  skips when OpenMP is unavailable.
- Keep tests deterministic across platforms (macOS/Windows/Linux all
  run in CI).

## Code Style

- R: exported functions follow upstream names (`trust`, `bobyqa`,
  `newuoa`, `steihaug`) to stay drop-in compatible; internal helpers use
  a leading dot (`.trustTryEval`, `.minqaCommonArgs`).
- C/C++: C entry points and types use `snake_case` with solver prefix
  (`trust_`, `minqa_`, `steihaug_`). Use `Rf_`-prefixed R API (the
  public header defines `R_NO_REMAP`).
- Keep argument semantics, defaults and return values identical to the
  upstream R packages ('trust', 'minqa') unless NEWS.md documents a
  deliberate deviation (e.g. BOBYQA's `RESCUE` step fix).
- ASCII only. No Unicode anywhere in the repo (CRAN requirement): use
  `--` for dashes, `->` for arrows, straight quotes, `...` for ellipses.

## Release / CRAN Notes

- Update `NEWS.md` for user-visible changes and `cran-comments.md` on
  submission.
- When porting or adding third-party code, add its source and notice to
  `inst/COPYRIGHTS` and credit authors in `DESCRIPTION` `Authors@R`.
- `.Rbuildignore` excludes `README.Rmd`, `_pkgdown.yml`, `pkgdown/`,
  `cran-comments.md`, `.github/`, `CLAUDE.md`, etc.; add new dev-only
  files there.
- The `RcppTrust.Rmd` vignette runs `microbenchmark` comparisons
  (conditional on Suggests being installed); it must only report
  timings, never assert on them.

## CRAN tests: never assert on timing

- Never put a wall-clock timing expectation in a test that can run on CRAN:
  no `expect_lt()`/`expect_true()` on `Sys.time()` differences,
  `system.time()`, `proc.time()`, `bench`/`microbenchmark` results, or any
  other elapsed-time threshold. CRAN check machines are shared and heavily
  loaded, so elapsed time there is random and a timing assertion fails at
  random. Example: nlmixr2plot `test-model-diagram.R` asserted a
  `modelGraph()` call took < 30 s; on CRAN it took 48.4 s and failed.
- Test the result, not the speed. If a speed regression guard is really
  needed, put it only in a test that starts with `skip_on_cran()` (or the
  package's own CRAN skip helper), and never in examples or vignettes.
- Tests that take more than a few seconds should also skip on CRAN.
