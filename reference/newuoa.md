# Unconstrained derivative-free optimization (NEWUOA)

Thread-safe C++ port of
[`minqa::newuoa()`](https://rdrr.io/pkg/minqa/man/newuoa.html):
minimizes a function of many variables by a trust region method that
forms quadratic models by interpolation, using M. J. D. Powell's NEWUOA
algorithm.

## Usage

``` r
newuoa(par, fn, control = list(), ...)
```

## Source

M. J. D. Powell's original Fortran 77 NEWUOA code, archived at
<https://github.com/libprima/prima/tree/main/fortran/original/newuoa>; R
interface following the minqa package.

## Arguments

- par:

  numeric vector of starting parameters.

- fn:

  function to be minimized; its first argument must be the parameter
  vector and it must return a scalar numeric value.

- control:

  a list of control settings: `npt`, `rhobeg`, `rhoend`, `iprint` and
  `maxfun`, as in
  [`minqa::newuoa()`](https://rdrr.io/pkg/minqa/man/newuoa.html).

- ...:

  further arguments passed to `fn`.

## Value

A list of class `c("newuoa", "minqa")` with components `par`, `fval`,
`feval`, `ierr` and `msg`.

## References

M. J. D. Powell (2006), "The NEWUOA software for unconstrained
optimization without derivatives", in Large-Scale Nonlinear
Optimization, Springer, 255-297.
[doi:10.1007/0-387-30065-1_16](https://doi.org/10.1007/0-387-30065-1_16)

## See also

\[bobyqa()\], \[minqa_c_api\]

## Examples

``` r
fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
newuoa(c(1, 2), fr)
#> parameter estimates: 1.00000116212176, 1.00000231682727 
#> objective: 1.35602908489591e-12 
#> number of function evaluations: 136 
```
