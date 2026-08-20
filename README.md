
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RcppTrust <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
<!-- badges: end -->

RcppTrust is a thread-safe C++ port of the trust-region optimizer in
Charles J. Geyer’s CRAN package
[`trust`](https://cran.r-project.org/package=trust), built for use in
[`nlmixr2est`](https://github.com/nlmixr2/nlmixr2est). It implements the
exact same algorithm – same Newton / easy-easy / hard-easy / hard-hard
trust-region subproblem, same termination criteria – and exposes it
three ways:

- **`trust()`**, a drop-in R replacement for `trust::trust()`, with the
  same arguments, defaults, and return value.
- **`trust_solve_c()`**, a thread-safe C entry point (a C objective
  function pointer plus a small options struct, no R API calls) safe to
  call from parallel C++ code such as an OpenMP loop.
- A header-only, positionally-indexed function-pointer table
  (`inst/include/RcppTrust.h`), so another package can call
  `trust_solve_c()` without linking against this package’s shared
  library – the same pattern `rxode2`/`n1qn1`/`lbfgsb3c` use across the
  nlmixr2 ecosystem.

See `vignette("RcppTrust")` for what’s the same as upstream `trust`,
what’s different, and a worked example of the C interface and the
registration pattern.

Note this package was generated with the help of AI (Claude/Gemini).

## Installation

You can install the development version of RcppTrust from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("nlmixr2/RcppTrust")
```

## Example

`trust()` is a straight substitute for `trust::trust()`:

``` r
library(RcppTrust)

# Rosenbrock's function, the example from ?trust::trust
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

## Authors

- Charles J. Geyer – original algorithm and R implementation (`trust`)
- Matthew Fidler – thread-safe C++ port (`RcppTrust`)
