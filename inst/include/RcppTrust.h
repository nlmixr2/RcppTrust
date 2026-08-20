#ifndef __RCPPTRUST_H__
#define __RCPPTRUST_H__

/* Header-only, positionally-indexed function-pointer registration for
 * RcppTrust's thread-safe C entry point, modeled exactly on the pattern
 * used by n1qn1c/lbfgsb3c/rxode2 in the nlmixr2 ecosystem (see
 * n1qn1c/inst/include/n1qn1c.h for the smallest reference example).
 *
 * A consuming package never links against RcppTrust's shared library.
 * Instead it:
 *   1. Adds `RcppTrust` to its DESCRIPTION's LinkingTo (for this header)
 *      and Imports (for the R-level .RcppTrustPtr() getter).
 *   2. In exactly one translation unit, includes this header wrapped in
 *      renaming #defines so the generated globals/init function get a
 *      package-unique name (avoids collisions if multiple consumers are
 *      loaded together), e.g.:
 *
 *        extern "C" {
 *        #define iniRcppTrustPtrs _mypkg_iniRcppTrustPtrs
 *          iniRcppTrust
 *        }
 *
 *      then declares/registers `_mypkg_iniRcppTrustPtrs` as a .Call entry
 *      in its own src/init.c, taking one SEXP argument.
 *   3. In R's .onLoad(), calls:
 *
 *        .Call(`_mypkg_iniRcppTrustPtrs`, RcppTrust:::.RcppTrustPtr(),
 *              PACKAGE = "mypkg")
 *
 *      which fills the function pointer(s) below from the positionally-
 *      indexed list RcppTrust::.RcppTrustPtr() returns.
 *
 * New slots are only ever appended at the end, never reordered or
 * removed, so upgrading RcppTrust never breaks a consumer built against
 * an older version (the documented nlmixr2est convention this mirrors). */

/* Avoid Rinternals.h's CamelCase-free macro remapping (length, error,
 * ...), which otherwise collides with C++ standard library member names
 * when this header is pulled in ahead of Rcpp/RcppArmadillo -- as
 * compileAttributes() does automatically for a header matching the
 * package name. Only Rf_-prefixed / CamelCase API is used below. */
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "trust_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

typedef int (*trust_solve_c_t)(int n, const double *parinit, trust_c_objfun_t objfun,
                                void *userdata, const trust_options_t *opts,
                                trust_result_t *result);
typedef void (*trust_result_free_t)(trust_result_t *res);

extern trust_solve_c_t trust_solve_c_ptr;
extern trust_result_free_t trust_result_free_ptr;

static inline SEXP iniRcppTrustPtrs0(SEXP p) {
  if (trust_solve_c_ptr == NULL) {
    trust_solve_c_ptr = (trust_solve_c_t)R_ExternalPtrAddrFn(VECTOR_ELT(p, 0));
    trust_result_free_ptr =
        (trust_result_free_t)R_ExternalPtrAddrFn(VECTOR_ELT(p, 1));
  }
  return R_NilValue;
}

#define iniRcppTrust                             \
  trust_solve_c_t trust_solve_c_ptr = NULL;       \
  trust_result_free_t trust_result_free_ptr = NULL; \
  SEXP iniRcppTrustPtrs(SEXP p) {                 \
    iniRcppTrustPtrs0(p);                         \
    return R_NilValue;                            \
  }

#if defined(__cplusplus)
}
#endif

#endif /* __RCPPTRUST_H__ */
