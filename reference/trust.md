# Non-Linear Optimization by a Trust Region Algorithm

Thread-safe C++ port of
[`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html) (Geyer).
Carries out a minimization or maximization of a function using a trust
region algorithm; see the original package's documentation for the
algorithm and argument semantics, which this function reproduces
exactly.

## Usage

``` r
trust(
  objfun,
  parinit,
  rinit,
  rmax,
  parscale,
  iterlim = 100,
  fterm = sqrt(.Machine$double.eps),
  mterm = sqrt(.Machine$double.eps),
  minimize = TRUE,
  blather = FALSE,
  ...
)
```

## Arguments

- objfun:

  an R function computing value, gradient, and Hessian; see
  [`trust::trust`](https://rdrr.io/pkg/trust/man/trust.html) for the
  exact contract.

- parinit:

  starting parameter values.

- rinit:

  starting trust region radius.

- rmax:

  maximum allowed trust region radius.

- parscale:

  optional parameter scale vector.

- iterlim:

  maximum number of iterations.

- fterm:

  function-value termination tolerance.

- mterm:

  model-decrease termination tolerance.

- minimize:

  if `TRUE` minimize, if `FALSE` maximize.

- blather:

  if `TRUE` return extra per-iteration information.

- ...:

  additional arguments passed to `objfun`.

## Value

A list with the same components as
[`trust::trust()`](https://rdrr.io/pkg/trust/man/trust.html).
