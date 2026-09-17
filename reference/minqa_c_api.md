# Thread-safe C interfaces

Every optimizer in this package is implemented as a C-callable core with
no R API calls and no global state, so each can run concurrently from
OpenMP/parallel C++ code (e.g. as an inner or outer optimizer in
`nlmixr2est`). The R functions \[trust()\], \[steihaug()\], \[bobyqa()\]
and \[newuoa()\] are thin wrappers around these same cores.

## Details

|            |                    |                                             |
|------------|--------------------|---------------------------------------------|
| **Solver** | **Header**         | **Entry point / free**                      |
| trust      | `trust_types.h`    | `trust_solve_c` / `trust_result_free`       |
| BOBYQA     | `minqa_types.h`    | `bobyqa_solve_c` / `minqa_result_free`      |
| NEWUOA     | `minqa_types.h`    | `newuoa_solve_c` / `minqa_result_free`      |
| Steihaug   | `steihaug_types.h` | `steihaug_solve_c` / `steihaug_result_free` |

A consuming package adds `RcppTrust` to `LinkingTo` and `Imports`,
includes `RcppTrust.h`, and resolves the pointers at load time from
\[.RcppTrustPtr()\] (slots, in order: `trust_solve_c`,
`trust_result_free`, `bobyqa_solve_c`, `newuoa_solve_c`,
`minqa_result_free`, `steihaug_solve_c`, `steihaug_result_free`). See
[`vignette("RcppTrust")`](../articles/RcppTrust.md) for a worked
example.
