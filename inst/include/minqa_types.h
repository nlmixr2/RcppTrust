#ifndef MINQA_TYPES_H
#define MINQA_TYPES_H

/* Plain C types for the thread-safe ports of M. J. D. Powell's
 * derivative-free quadratic-model optimizers BOBYQA and NEWUOA,
 * as distributed (with R-specific edits) in the CRAN package 'minqa'.
 * No R API and no C++ types appear here so this header can be included
 * from a consumer's plain C code as well. */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Objective function contract for the thread-safe minqa solvers.
 *
 * n        : number of parameters
 * x        : parameter vector, length n (input, never non-finite)
 * f        : output, the objective value at x. A non-finite value is
 *            replaced by DBL_MAX exactly as minqa's calfun() does.
 * userdata : opaque pointer passed through unchanged
 *
 * Return value: 0 on success, <0 to abort the optimization (the solver
 * returns immediately with result->ierr == MINQA_ERR_OBJFUN). */
typedef int (*minqa_c_objfun_t)(int n, const double *x, double *f,
                                void *userdata);

/* Optional sink for iprint output. Called with one NUL-terminated,
 * newline-terminated line at a time. NULL disables all printing (the
 * default for the thread-safe path, where Rprintf() must not be used). */
typedef void (*minqa_print_t)(const char *line, void *print_userdata);

typedef struct {
  int npt;       /* number of interpolation points */
  double rhobeg; /* initial trust region radius, > 0 */
  double rhoend; /* final trust region radius, > 0 (normally <= rhobeg) */
  int iprint;    /* 0..3 as in minqa (>3: print every iprint-th eval) */
  int maxfun;    /* maximum number of objective evaluations */
  minqa_print_t print;  /* NULL = no printing */
  void *print_userdata;
} minqa_options_t;

/* Initializer mirroring minqa's commonArgs() defaults for a problem with
 * n parameters starting at par (rhobeg = min(0.95, 0.2 * max(abs(par))),
 * rhoend = 1e-6 * rhobeg -- both 0, hence invalid, when every par is 0,
 * exactly as minqa rejects that start -- npt = max(n + 2, min(min(n + 2, 2n),
 * (n + 1)(n + 2) / 2)), iprint = 0, maxfun = 10000). */
static inline minqa_options_t minqa_options_default(int n, const double *par) {
  minqa_options_t opts;
  double mx = 0.0;
  int i, npt, nptmax;
  for (i = 0; i < n; i++) {
    double a = par[i] < 0 ? -par[i] : par[i];
    if (i == 0 || a > mx) mx = a;
  }
  opts.rhobeg = 0.2 * mx;
  if (opts.rhobeg > 0.95) opts.rhobeg = 0.95;
  opts.rhoend = 1.0e-6 * opts.rhobeg;
  npt = n + 2 < 2 * n ? n + 2 : 2 * n;
  nptmax = ((n + 1) * (n + 2)) / 2;
  if (npt > nptmax) npt = nptmax;
  if (npt < n + 2) npt = n + 2;
  opts.npt = npt;
  opts.iprint = 0;
  opts.maxfun = 10000;
  opts.print = NULL;
  opts.print_userdata = NULL;
  return opts;
}

/* result->ierr codes. Non-negative codes are exactly the codes minqa's
 * R wrappers report (after their remapping of the raw Fortran IERR):
 *   0 normal exit
 *   1 maximum number of function evaluations exceeded   (Fortran 390)
 *   2 NPT is not in the required interval               (Fortran 10)
 *   3 a trust region step failed to reduce q            (430/3701)
 *   4 a box constraint range is too small (< 2*rhobeg)  (Fortran 20)
 *   5 too much cancellation in a denominator (BOBYQA only, Fortran 320)
 * Negative codes are failures specific to this C interface. */
#define MINQA_OK 0
#define MINQA_ERR_MAXFUN 1
#define MINQA_ERR_NPT 2
#define MINQA_ERR_NO_REDUCTION 3
#define MINQA_ERR_RANGE 4
#define MINQA_ERR_CANCELLATION 5
#define MINQA_ERR_OBJFUN -1    /* objfun returned < 0 */
#define MINQA_ERR_NONFINITE_X -2 /* calfun received non-finite x (minqa stops) */
#define MINQA_ERR_INVALID -3 /* n < 0, NULL par/objfun, or rhobeg/rhoend
                             * not both > 0. rhoend > rhobeg is allowed, as
                             * minqa's bobyqa() can itself produce it when it
                             * shrinks rhobeg for tight bounds. */
#define MINQA_ERR_NOMEM -4
#define MINQA_ERR_INTERNAL -99

typedef struct {
  int n;
  int ierr;        /* see codes above */
  int raw_ierr;    /* the untranslated Fortran IERR (0, 10, 20, 320, ...) */
  int feval;       /* number of objective evaluations made by the solver */
  double fval;     /* objective at par, re-evaluated after the solver
                    * returns (not counted in feval, not DBL_MAX-clamped),
                    * as minqa's rval() does */
  double *par;     /* length n, owned by this struct (malloc) */
} minqa_result_t;

/* Zero-initializes a minqa_result_t (does not allocate). */
void minqa_result_zero(minqa_result_t *res);

/* Frees the buffers owned by *res (safe on a zeroed result) and zeroes it. */
void minqa_result_free(minqa_result_t *res);

/* Thread-safe solvers: no R API, no global/static mutable state; each call
 * allocates its own workspace. Return value == result->ierr
 * (result == NULL returns MINQA_ERR_INVALID without writing anything).
 *
 * opts may be NULL, meaning minqa_options_default(n, par). For bobyqa,
 * lower/upper may be NULL, meaning -Inf/+Inf for every parameter.
 *
 * Deviation from minqa 1.2.8: minqa's rescue.f calls CALFUN(N,X,IPRINT)
 * with X undeclared in RESCUE (an uninitialized scalar), so BOBYQA's
 * RESCUE step evaluates the objective at garbage memory. This port uses
 * Powell's original CALFUN(N,W,F), i.e. the intended point W(1..N). All
 * runs that never enter RESCUE match minqa bit for bit. */
int bobyqa_solve_c(int n, const double *par, const double *lower,
                   const double *upper, minqa_c_objfun_t objfun,
                   void *userdata, const minqa_options_t *opts,
                   minqa_result_t *result);

int newuoa_solve_c(int n, const double *par, minqa_c_objfun_t objfun,
                   void *userdata, const minqa_options_t *opts,
                   minqa_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* MINQA_TYPES_H */
