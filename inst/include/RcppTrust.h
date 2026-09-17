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
#include "minqa_types.h"
#include "steihaug_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

typedef int (*trust_solve_c_t)(int n, const double *parinit, trust_c_objfun_t objfun,
                                void *userdata, const trust_options_t *opts,
                                trust_result_t *result);
typedef void (*trust_result_free_t)(trust_result_t *res);
typedef int (*bobyqa_solve_c_t)(int n, const double *par, const double *lower,
                                const double *upper, minqa_c_objfun_t objfun,
                                void *userdata, const minqa_options_t *opts,
                                minqa_result_t *result);
typedef int (*newuoa_solve_c_t)(int n, const double *par,
                                minqa_c_objfun_t objfun, void *userdata,
                                const minqa_options_t *opts,
                                minqa_result_t *result);
typedef void (*minqa_result_free_t)(minqa_result_t *res);
typedef int (*steihaug_solve_c_t)(int n, const double *parinit,
                                  steihaug_c_objfun_t objfun,
                                  steihaug_c_hessvec_t hessvec_or_null,
                                  void *userdata,
                                  const steihaug_options_t *opts,
                                  steihaug_result_t *result);
typedef void (*steihaug_result_free_t)(steihaug_result_t *res);

extern trust_solve_c_t trust_solve_c_ptr;
extern trust_result_free_t trust_result_free_ptr;
extern bobyqa_solve_c_t bobyqa_solve_c_ptr;
extern newuoa_solve_c_t newuoa_solve_c_ptr;
extern minqa_result_free_t minqa_result_free_ptr;
extern steihaug_solve_c_t steihaug_solve_c_ptr;
extern steihaug_result_free_t steihaug_result_free_ptr;

/* Slots beyond the list's length (an older RcppTrust) are left NULL, so a
 * consumer should check e.g. `steihaug_solve_c_ptr != NULL` before use. */
static inline DL_FUNC iniRcppTrustSlot0(SEXP p, R_xlen_t i) {
  return Rf_xlength(p) > i ? R_ExternalPtrAddrFn(VECTOR_ELT(p, i)) : NULL;
}

static inline SEXP iniRcppTrustPtrs0(SEXP p) {
  if (trust_solve_c_ptr == NULL) {
    trust_solve_c_ptr = (trust_solve_c_t)iniRcppTrustSlot0(p, 0);
    trust_result_free_ptr = (trust_result_free_t)iniRcppTrustSlot0(p, 1);
    bobyqa_solve_c_ptr = (bobyqa_solve_c_t)iniRcppTrustSlot0(p, 2);
    newuoa_solve_c_ptr = (newuoa_solve_c_t)iniRcppTrustSlot0(p, 3);
    minqa_result_free_ptr = (minqa_result_free_t)iniRcppTrustSlot0(p, 4);
    steihaug_solve_c_ptr = (steihaug_solve_c_t)iniRcppTrustSlot0(p, 5);
    steihaug_result_free_ptr = (steihaug_result_free_t)iniRcppTrustSlot0(p, 6);
  }
  return R_NilValue;
}

#define iniRcppTrust                                    \
  trust_solve_c_t trust_solve_c_ptr = NULL;              \
  trust_result_free_t trust_result_free_ptr = NULL;      \
  bobyqa_solve_c_t bobyqa_solve_c_ptr = NULL;            \
  newuoa_solve_c_t newuoa_solve_c_ptr = NULL;            \
  minqa_result_free_t minqa_result_free_ptr = NULL;      \
  steihaug_solve_c_t steihaug_solve_c_ptr = NULL;        \
  steihaug_result_free_t steihaug_result_free_ptr = NULL; \
  SEXP iniRcppTrustPtrs(SEXP p) {                        \
    iniRcppTrustPtrs0(p);                                \
    return R_NilValue;                                   \
  }

#if defined(__cplusplus)
}
#endif

#endif /* __RCPPTRUST_H__ */
