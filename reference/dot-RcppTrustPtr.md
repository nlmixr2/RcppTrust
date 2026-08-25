# Positionally-indexed function pointer table for downstream packages

Returns the list of native function pointers a consuming package (e.g.
nlmixr2est) resolves via the header-only registration pattern in
`inst/include/RcppTrust.h` – the same pattern rxode2/n1qn1c/ lbfgsb3c
use, so no linker dependency on this package's shared library is
required. Not intended to be called directly by end users.

## Usage

``` r
.RcppTrustPtr()
```

## Value

An unclassed list of external pointers.
