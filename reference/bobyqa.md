# Nonlinear optimization with box constraints (BOBYQA)

Thread-safe C++ port of
[`minqa::bobyqa()`](https://rdrr.io/pkg/minqa/man/bobyqa.html):
minimizes a function of many variables subject to box constraints by a
trust region method that forms quadratic models by interpolation, using
M. J. D. Powell's BOBYQA algorithm. Arguments, defaults, and return
value are those of
[`minqa::bobyqa()`](https://rdrr.io/pkg/minqa/man/bobyqa.html); see its
documentation for details.

## Usage

``` r
bobyqa(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
```

## Source

M. J. D. Powell's original Fortran 77 BOBYQA code, archived at
<https://github.com/libprima/prima/tree/main/fortran/original/bobyqa>; R
interface following the minqa package.

## Arguments

- par:

  numeric vector of starting parameters.

- fn:

  function to be minimized; its first argument must be the parameter
  vector and it must return a scalar numeric value.

- lower, upper:

  numeric vectors of lower and upper bounds (recycled when of length 1).

- control:

  a list of control settings: `npt`, `rhobeg`, `rhoend`, `iprint`,
  `maxfun`, `obstop` and `force.start`, as in
  [`minqa::bobyqa()`](https://rdrr.io/pkg/minqa/man/bobyqa.html).

- ...:

  further arguments passed to `fn`.

## Value

A list of class `c("bobyqa", "minqa")` with components `par`, `fval`,
`feval`, `ierr` and `msg`.

## References

M. J. D. Powell (2009), "The BOBYQA algorithm for bound constrained
optimization without derivatives", Report No. DAMTP 2009/NA06, Centre
for Mathematical Sciences, University of Cambridge.

## See also

\[newuoa()\], and \[minqa_c_api\] for calling the solvers from parallel
C/C++ code.

## Examples

``` r
fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4))
#> parameter estimates: 0.999999968901681, 0.999999928305543 
#> objective: 9.98796325533312e-15 
#> number of function evaluations: 341 
```
