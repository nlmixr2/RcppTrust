# RcppTrust: A Thread-Safe C++ Trust Region Optimizer

`RcppTrust` is a C++ port of the trust-region optimizer in Charles J.
Geyer’s CRAN package `trust`, written for use inside `nlmixr2est`. This
vignette does not re-derive the algorithm – see Geyer’s original
paper/vignette, or Nocedal and Wright (1999, Chapter 4) and Fletcher
(1987, Section 5.1), for that. Instead it covers what a user of either
package actually needs to know: what carries over unchanged, what’s new,
and how to call the new thread-safe C interface from your own package.

## Similarities: what carries over unchanged

[`RcppTrust::trust()`](../reference/trust.md) is a drop-in replacement
for [`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html). It has
the same arguments, the same defaults, the same return value, and
implements the *exact same algorithm*: the same Newton / easy-easy /
hard-easy / hard-hard case split for the trust-region subproblem, the
same accept/reject and radius-adjustment rule, and the same termination
criteria.

``` r

library(RcppTrust)

# Rosenbrock's function, exactly the example from ?trust::trust
objfun <- function(x) {
  f <- expression(100 * (x2 - x1^2)^2 + (1 - x1)^2)
  g1 <- D(f, "x1"); g2 <- D(f, "x2")
  h11 <- D(g1, "x1"); h12 <- D(g1, "x2"); h22 <- D(g2, "x2")
  x1 <- x[1]; x2 <- x[2]
  list(
    value = eval(f), gradient = c(eval(g1), eval(g2)),
    hessian = rbind(c(eval(h11), eval(h12)), c(eval(h12), eval(h22)))
  )
}

out <- trust(objfun, c(3, 1), 1, 5)
out[c("value", "argument", "converged", "iterations")]
#> $value
#> [1] 5.165437e-15
#> 
#> $argument
#> [1] 1 1
#> 
#> $converged
#> [1] TRUE
#> 
#> $iterations
#> [1] 21
```

Everything you already know about
[`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html) – the
`objfun` contract (return a list with `value`, `gradient`, `hessian`, or
`list(value = Inf)`/`list(value = -Inf)` to signal an infeasible point),
the `rinit`/`rmax`/`parscale`/`iterlim`/`fterm`/`mterm`/
`minimize`/`blather` arguments, the shape of the returned list,
`blather = TRUE`’s extra `argpath`/`steptype`/`rho`/… components –
applies unchanged. All 10 of upstream’s own `tests/*.R` scripts are
ported into this package’s test suite and pass against the C++ core.

## Differences

|  | `trust` | `RcppTrust` |
|----|----|----|
| Implementation | pure R | C++ (RcppArmadillo for the linear algebra) |
| Root-finder for the trust-region subproblem | [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) | `boost::math::tools::toms748_solve()` |
| Callable from R | yes | yes, identically ([`trust()`](../reference/trust.md)) |
| Callable from thread-safe C/C++ | no | yes (`trust_solve_c()`, see below) |
| Usable from another package without linking against its shared library | no C API at all | yes, via a header-only function-pointer table |
| Authors | Charles J. Geyer | Charles J. Geyer (algorithm), Matthew Fidler (C++ port) |

A few of these are worth expanding on:

- **Numeric agreement, not bit-identical output.** Swapping
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) for
  `toms748_solve()` (both are derivative-free, bracketing root-finders,
  just different implementations) means results agree with upstream to
  root-finder tolerance (`.Machine$double.eps^0.25`, upstream’s own
  default) rather than to the last bit – typically 1e-4 to 1e-8 in the
  returned `argument`, accumulating slightly over iterations.
  Convergence, iteration count, and the sequence of `steptype`s match
  exactly.

- **The hard-hard case can mirror-flip a sign.** When the trust-region
  subproblem hits an exact eigenvalue tie (the “hard-hard” case), the
  solution involves an eigenvector that’s only defined up to sign.
  RcppTrust uses RcppArmadillo/LAPACK’s `eig_sym()` rather than R’s
  [`eigen()`](https://rdrr.io/r/base/eigen.html), which can pick the
  opposite sign than upstream did. This produces a mirror-image (equally
  valid) solution with the *same* objective value, convergence, and
  iteration count – just a flipped sign on the affected component. This
  is a property of eigendecomposition at exact ties, not a bug in the
  port.

- **Thread safety is new.**
  [`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html) calls
  back into R on every iteration and so can only ever run on R’s main
  thread. Nothing about [`RcppTrust::trust()`](../reference/trust.md)
  changes that – it *also* calls back into R for the objective function,
  and is not more thread-safe than upstream. What’s new is a second,
  separate entry point, `trust_solve_c()`, that never touches R at all
  and is safe to call from parallel C++ code (verified with an OpenMP
  stress test in this package’s own test suite).

- **Iteration counts can diverge, occasionally by a lot, on
  ill-conditioned problems.** Both implementations always converge to
  the same solution, but the *path* they take to get there isn’t
  guaranteed to match: a case split (typically easy-easy vs. hard-easy)
  can tip the other way from a tiny eigendecomposition difference
  between R’s [`eigen()`](https://rdrr.io/r/base/eigen.html) and
  LAPACK’s `eig_sym()`, changing which steps get accepted and how the
  trust region radius evolves for the rest of the run. The
  [Speed](#speed) section below has a concrete example converging to the
  same point in 16 iterations (RcppTrust) vs. 28 (upstream).

## Speed

Both implementations solve the same trust-region subproblem, so this is
really comparing three things at once: R vs. C++ overhead per iteration,
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html)
vs. `boost::math::tools::toms748_solve()`, and (per the note above)
however many iterations each happens to take to converge on a given
problem. [`RcppTrust::trust()`](../reference/trust.md) still calls back
into R for the objective function on every iteration – exactly like
upstream – so none of this is about avoiding R call overhead; it’s the
trust-region bookkeeping itself (the eigendecomposition and root-find)
that’s faster in C++, and that saving grows with the number of
parameters.

``` r

library(microbenchmark)

# same Rosenbrock objfun as above: 2 parameters, ~20 iterations
mb_small <- microbenchmark(
  trust = trust::trust(objfun, c(3, 1), 1, 5),
  RcppTrust = RcppTrust::trust(objfun, c(3, 1), 1, 5),
  times = 50
)
print(mb_small)
#> Unit: microseconds
#>       expr      min       lq      mean    median       uq      max neval
#>      trust 2135.970 2182.005 2389.2443 2226.1225 2316.006 4647.711    50
#>  RcppTrust  817.987  853.464  880.8299  875.0095  912.084  978.438    50
```

A larger problem – the restricted-domain log-barrier objective from
[`?trust::trust`](https://rdrr.io/pkg/trust/man/trust.html)’s second
example, at `d = 30` parameters instead of 5 – makes the per-iteration
bookkeeping cost (rather than the objective function itself) a bigger
share of the total, and is the case behind the 28-vs-16-iteration
example mentioned above:

``` r

d <- 30
mu <- seq_len(d)
barrier_objfun <- function(x) {
  normxsq <- sum(x^2)
  omnormxsq <- 1 - normxsq
  if (normxsq >= 1) return(list(value = Inf))
  f <- sum(x * mu) - log(omnormxsq)
  g <- mu + 2 * x / omnormxsq
  B <- 4 * outer(x, x) / omnormxsq^2 + 2 * diag(d) / omnormxsq
  list(value = f, gradient = g, hessian = B)
}

r1 <- trust::trust(barrier_objfun, rep(0, d), 1, 100)
r2 <- RcppTrust::trust(barrier_objfun, rep(0, d), 1, 100)
# same solution, different number of steps to get there (see above)
c(trust_iterations = r1$iterations, RcppTrust_iterations = r2$iterations)
#>     trust_iterations RcppTrust_iterations 
#>                   32                   16
max(abs(r1$argument - r2$argument))
#> [1] 3.989087e-12

mb_large <- microbenchmark(
  trust = trust::trust(barrier_objfun, rep(0, d), 1, 100),
  RcppTrust = RcppTrust::trust(barrier_objfun, rep(0, d), 1, 100),
  times = 30
)
print(mb_large)
#> Unit: microseconds
#>       expr      min       lq     mean   median       uq       max neval
#>      trust 5799.622 5851.359 6649.812 5931.348 7443.403 11499.358    30
#>  RcppTrust  914.498 1002.923 1088.400 1059.224 1183.150  1331.357    30
```

On this machine, [`RcppTrust::trust()`](../reference/trust.md) comes out
roughly 2x faster on the 2-parameter problem and roughly 4-5x faster on
the 30-parameter one, despite taking *more* per-iteration R/C++ round
trips than upstream (each iteration’s objective-function call
additionally goes through the internal `.trustTryEval()` wrapper
described below, so that objfun errors are still caught exactly the way
upstream’s own [`try()`](https://rdrr.io/r/base/try.html) catches them)
– the win comes entirely from the trust-region bookkeeping itself, and
should be expected to grow with the number of parameters. It says
nothing about the thread-safe C path (`trust_solve_c()`), which
additionally removes the R round trip altogether; see below.

## When to use which interface

- **Writing an R script, or need a drop-in replacement for
  [`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html)?** Use
  [`RcppTrust::trust()`](../reference/trust.md). It’s a straight
  substitute.
- **Calling from your own C/C++ code, especially from inside a parallel
  loop (OpenMP, or similar)?** Use the C interface, `trust_solve_c()`,
  described below. This is the intended use in `nlmixr2est`,
  e.g. fitting each subject’s inner problem on its own thread.
- **Writing a package that wants to call `trust_solve_c()` without
  taking on a hard link-time dependency on RcppTrust’s shared library?**
  Use the header-only registration mechanism, also described below –
  it’s how the C interface is actually meant to be consumed by another
  package.

## Using the C interface

The thread-safe entry point is declared in `inst/include/trust_types.h`
(a plain C header with no R and no C++ dependency, so it’s safe to
include from plain C code too):

``` c
typedef int (*trust_c_objfun_t)(int n, const double *par,
    double *value, double *gradient, double *hessian, void *userdata);
```

Your objective function fills `*value` (and `gradient`/`hessian`, each
length `n`/`n*n`, unless the point is infeasible) and returns:

- `0` – feasible: `*value`, `gradient`, `hessian` are all filled and
  finite.
- `1` – infeasible/out of the objective’s domain: only `*value` is set,
  to `+INFINITY` when minimizing or `-INFINITY` when maximizing
  (`gradient`/`hessian` are left untouched). This is the C-level
  equivalent of upstream’s `list(value = Inf)` convention.
- negative – a hard error (the evaluation could not be completed at
  all).

Options are a plain struct, built from `trust_options_default()` (which
fills in [`trust()`](../reference/trust.md)’s R-level defaults) and then
adjusted as needed:

``` c
trust_options_t opts = trust_options_default(/* rinit = */ 1.0, /* rmax = */ 5.0);
opts.iterlim = 200;      // default 100
opts.minimize = 0;       // maximize instead
opts.blather = 1;        // collect per-iteration arrays, like blather = TRUE
opts.has_parscale = 1;
opts.parscale = my_parscale;  // a `const double *` of length n you own
```

The result is written into a caller-supplied `trust_result_t`, which
owns its buffers (`argument`, `gradient`, `hessian`, and – if
`opts.blather` – the `argpath`/`argtry`/`steptype`/`rho`/… arrays, laid
out exactly like the R-level `blather = TRUE` output) until you free it:

``` c
trust_result_t res;
int rc = trust_solve_c(n, parinit, my_objfun, my_userdata, &opts, &res);
// rc == 0 on success; res.converged, res.iterations, res.argument, ...
trust_result_free(&res);
```

`trust_solve_c()` and `trust_result_free()` take no locks and touch no
shared, mutable state – every argument is either an input or owned
exclusively by the caller’s `trust_result_t` – so many threads can each
be running their own `trust_solve_c()` call at the same time, on
independent problems, with no coordination needed. That’s the piece
[`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html) cannot
offer at all.

### Calling it from your own package: the registration pattern

`RcppTrust`’s shared library is never linked against directly. Instead –
following the same header-only, positionally-indexed function-pointer
pattern already used across the nlmixr2 ecosystem (`rxode2`, `n1qn1`,
`lbfgsb3c`) – a consumer resolves `trust_solve_c`/`trust_result_free` as
function pointers once, at load time. The full wiring is three pieces,
which the example below reproduces in a single self-contained file for
demonstration.

1.  Add `RcppTrust` to your `DESCRIPTION`’s `LinkingTo` (for the header)
    and `Imports` (for the R-level
    [`.RcppTrustPtr()`](../reference/dot-RcppTrustPtr.md) getter).

2.  In one translation unit, include the registration header with a
    package-unique renaming `#define` so the generated init function
    doesn’t collide with any other consumer’s:

    ``` cpp
    extern "C" {
    #define iniRcppTrustPtrs _mypkg_iniRcppTrustPtrs
    #include <RcppTrust.h>
    iniRcppTrust
    }
    ```

3.  In your package’s `.onLoad()`, resolve the pointers once:

    ``` r

    .onLoad <- function(libname, pkgname) {
      .Call(`_mypkg_iniRcppTrustPtrs`, RcppTrust:::.RcppTrustPtr(), PACKAGE = "mypkg")
    }
    ```

    After that, `trust_solve_c_ptr` and `trust_result_free_ptr` are live
    function pointers with the same signatures as
    `trust_solve_c()`/`trust_result_free()`, usable anywhere in your
    package’s C++ – including inside an OpenMP loop fitting many
    subjects’ problems in parallel:

    ``` cpp
    #pragma omp parallel for
    for (int i = 0; i < nSubjects; i++) {
      trust_options_t opts = trust_options_default(1.0, 5.0);
      trust_result_t res;
      trust_solve_c_ptr(n, parinit[i], subject_objfun, &subjectData[i], &opts, &res);
      // ... use res.argument, res.value, res.converged ...
      trust_result_free_ptr(&res);
    }
    ```

The example below is the whole pattern collapsed into one file and run
live, via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html),
exactly as it would work split across a real package’s `src/init.c`,
`R/zzz.R`, and wherever the fit happens:

``` r

cpp_code <- '
// [[Rcpp::depends(RcppTrust)]]
#include <Rcpp.h>

extern "C" {
#define iniRcppTrustPtrs _vignette_iniRcppTrustPtrs
#include <RcppTrust.h>
iniRcppTrust
}

// A thread-safe C objective function: Rosenbrock again, this time
// filling value/gradient/hessian directly instead of returning a list.
extern "C" int rosenbrock_c(int n, const double *par, double *value,
                             double *gradient, double *hessian, void *ud) {
  double x1 = par[0], x2 = par[1];
  double t = x2 - x1 * x1;
  *value = 100.0 * t * t + (1.0 - x1) * (1.0 - x1);
  gradient[0] = -400.0 * x1 * t - 2.0 * (1.0 - x1);
  gradient[1] = 200.0 * t;
  hessian[0] = -400.0 * x2 + 1200.0 * x1 * x1 + 2.0;
  hessian[1] = hessian[2] = -400.0 * x1;
  hessian[3] = 200.0;
  return 0;
}

// [[Rcpp::export]]
Rcpp::List fit_rosenbrock(SEXP ptrTable, Rcpp::NumericVector parinit) {
  // normally done once, in .onLoad() -- see above
  _vignette_iniRcppTrustPtrs(ptrTable);

  trust_options_t opts = trust_options_default(1.0, 5.0);
  trust_result_t res;
  trust_solve_c_ptr(parinit.size(), parinit.begin(), rosenbrock_c, nullptr, &opts, &res);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["argument"] = Rcpp::NumericVector(res.argument, res.argument + res.n),
    Rcpp::_["value"] = res.value,
    Rcpp::_["converged"] = res.converged != 0,
    Rcpp::_["iterations"] = res.iterations);
  trust_result_free_ptr(&res);
  return out;
}
'
Rcpp::sourceCpp(code = cpp_code)

fit_rosenbrock(RcppTrust:::.RcppTrustPtr(), c(3, 1))
#> $argument
#> [1] 1 1
#> 
#> $value
#> [1] 5.165437e-15
#> 
#> $converged
#> [1] TRUE
#> 
#> $iterations
#> [1] 21
```

Note that this compiled file never links against `RcppTrust`’s shared
library at all – `LinkingTo`/`Rcpp::depends()` only adds its
`inst/include` directory to the compiler’s include path. Every `trust_*`
symbol used above except the two resolved pointers (`trust_solve_c_ptr`,
`trust_result_free_ptr`) is a type or macro, not a function call, which
is exactly the point: nothing here creates an ABI dependency on a
specific build of RcppTrust, only on the stable, append-only pointer
table it publishes.
