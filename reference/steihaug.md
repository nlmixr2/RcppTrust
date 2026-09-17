# Steihaug truncated-CG trust-region Newton minimization

Thread-safe C++ port of the trust-region Newton minimizer with the
Steihaug truncated conjugate-gradient subproblem from the Rust crate
[basin](https://github.com/jolars/basin) (`TrustRegion` with its default
`Steihaug` subproblem, driven by basin's executor). Each iteration
approximately minimizes the quadratic model \\m(p) = f + g^T p +
\frac{1}{2} p^T B p\\ over \\\\p\\ \le \Delta\\ by conjugate gradients,
stopping at the boundary on negative curvature. Unlike \[trust()\] it
never factorizes the Hessian, and with `hessvec` it never forms it at
all (basin's matrix-free mode).

## Usage

``` r
steihaug(par, fn, gr, hess = NULL, hessvec = NULL, control = list(), ...)
```

## Arguments

- par:

  numeric vector of starting parameters.

- fn:

  objective function of the parameter vector returning a scalar.

- gr:

  gradient function returning a numeric vector of length `length(par)`.

- hess:

  Hessian function returning a symmetric `length(par)` by `length(par)`
  matrix. May be `NULL` when `hessvec` is supplied.

- hessvec:

  optional Hessian-vector product `function(x, v)` returning \\H(x) v\\.
  When supplied, the matrix-free mode is used and `hess` is never
  called.

- control:

  a list of control settings (defaults in brackets; a `NULL` tolerance
  or budget is disabled):

  rinit

  :   initial trust region radius \[1\]

  rmax

  :   maximum trust region radius \[100\]

  eta

  :   step acceptance threshold on the reduction ratio, in \\\[0, 1/4)\\
      \[0.125\]

  maxInner

  :   radius reductions (re-solves) per outer iteration \[10\]

  cgMaxit

  :   CG iteration cap per subproblem; 0 means `length(par)` \[0\]

  kappa, theta

  :   CG forcing rule: stop when \\\\r\\ \< \min(\kappa,
      \\g\\^\theta)\\g\\\\ \[0.5, 0.5\]

  maxit

  :   maximum outer iterations \[100\]

  maxCostEvals, maxGradEvals

  :   evaluation budgets \[NULL\]

  targetCost

  :   stop once the objective is at or below this \[NULL\]

  gradTol

  :   absolute gradient tolerance on \\\\g\\\_2\\ \[1e-8\]

  relGradTol

  :   gradient tolerance relative to the initial gradient norm \[NULL\]

  stepTol, relStepTol

  :   absolute / relative step tolerance \[NULL\]

  costTol, relCostTol

  :   absolute / relative objective change tolerance \[NULL\]

  trace

  :   if `TRUE`, return a per-attempt trace \[FALSE\]

- ...:

  further arguments passed to `fn`, `gr`, `hess` and `hessvec`.

## Value

A list with components `par`, `value`, `gradient`, `iterations`,
`converged`, `status` (integer), `message` (the basin termination
reason, e.g. `"GradientTolerance"`, `"MaxIter"`), `radius` (final trust
radius), `counts` (evaluations of `fn`, `gr`, `hess`, `hessvec`), and,
when `control$trace` is `TRUE`, `trace`, a data frame with one row per
subproblem attempt.

## Details

Evaluation pattern (identical to basin): `fn` and `gr` once at the
start; `hess` once per outer iteration (exact mode) or `hessvec` once
per CG product (matrix-free mode); `fn` at each trial point; `gr` after
each accepted step. `fn` may return `Inf` (or `NaN`) at a trial point to
reject it.

## References

Steihaug, T. (1983). The conjugate gradient method and trust regions in
large scale optimization. SIAM Journal on Numerical Analysis, 20(3),
626-637.

Nocedal, J. and Wright, S. J. (2006). Numerical Optimization (2nd ed.),
Algorithms 4.1 and 7.2. Springer.

## See also

\[trust()\], \[bobyqa()\], \[newuoa()\]

## Examples

``` r
fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
grr <- function(x) c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
                     200 * (x[2] - x[1]^2))
hr <- function(x) matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
                           -400 * x[1], 200), 2, 2)
steihaug(c(-1.2, 1), fr, grr, hr)
#> $par
#> [1] 1 1
#> 
#> $value
#> [1] 8.628062e-25
#> 
#> $gradient
#> [1]  3.688916e-11 -1.851852e-11
#> 
#> $iterations
#> [1] 27
#> 
#> $converged
#> [1] TRUE
#> 
#> $status
#> [1] 4
#> 
#> $message
#> [1] "GradientTolerance"
#> 
#> $radius
#> [1] 0.05108575
#> 
#> $counts
#>      fn      gr    hess hessvec 
#>      31      28      27       0 
#> 

# matrix-free: only Hessian-vector products
hv <- function(x, v) drop(hr(x) %*% v)
steihaug(c(-1.2, 1), fr, grr, hessvec = hv)$counts
#>      fn      gr    hess hessvec 
#>      31      28       0      76 
```
