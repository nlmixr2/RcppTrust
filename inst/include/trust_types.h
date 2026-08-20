#ifndef TRUST_TYPES_H
#define TRUST_TYPES_H

#include <stddef.h>

/* Plain C types shared by the thread-safe core, the registration header,
 * and the R interface layer. No R API and no C++ types appear here so
 * this file can be included from a consumer's plain C code as well. */

#ifdef __cplusplus
extern "C" {
#endif

/* Objective function contract for the thread-safe (C) path.
 *
 * n         : number of parameters
 * par       : parameter vector, length n (input)
 * value     : output, the objective value at par
 * gradient  : output buffer, length n (filled unless infeasible)
 * hessian   : output buffer, length n*n, row-major (filled unless infeasible)
 * userdata  : opaque pointer passed through unchanged
 *
 * Return value:
 *   0  feasible: *value, gradient[], hessian[] all filled and finite
 *   1  infeasible/out-of-domain: only *value is set, to +INFINITY when
 *      minimizing or -INFINITY when maximizing; gradient/hessian untouched
 *  <0  hard error (evaluation could not be completed at all)
 */
typedef int (*trust_c_objfun_t)(int n, const double *par,
                                 double *value, double *gradient,
                                 double *hessian, void *userdata);

typedef struct {
  double rinit;
  double rmax;
  int has_parscale;      /* 0/1 */
  const double *parscale; /* length n when has_parscale, else ignored */
  int iterlim;
  double fterm;
  double mterm;
  int minimize; /* 1 = minimize, 0 = maximize */
  int blather;  /* 1 = collect per-iteration blather arrays */
} trust_options_t;

/* Convenience initializer matching trust()'s R defaults, except rinit/rmax
 * which have no default upstream and must always be supplied. */
static inline trust_options_t trust_options_default(double rinit,
                                                      double rmax) {
  trust_options_t opts;
  opts.rinit = rinit;
  opts.rmax = rmax;
  opts.has_parscale = 0;
  opts.parscale = NULL;
  opts.iterlim = 100;
  opts.fterm = 1.4901161193847656e-08;  /* sqrt(.Machine$double.eps) */
  opts.mterm = 1.4901161193847656e-08;
  opts.minimize = 1;
  opts.blather = 0;
  return opts;
}

typedef struct {
  int n;
  /* 0 = no error. Positive: an objfun-raised error was caught gracefully
   * (mirrors upstream's try()-based recovery) -- 1 = on the very first
   * evaluation (no `argument`/`converged` beyond parinit is meaningful),
   * 2 = mid-loop on a trial point (`argument` holds that trial point),
   * 3 = only on the final post-loop re-evaluation (`argument`/`converged`
   * from the completed loop are still valid). Negative: a hard,
   * *uncaught* failure -- -2 objfun's return value failed validation
   * (check.objfun.output()-equivalent), -3 parinit itself was
   * infeasible, -4 out of memory, -99 an unexpected internal error. */
  int error;
  int converged;  /* 0/1 */
  int iterations;

  double value;
  double *gradient;  /* length n, owned by this struct */
  double *hessian;   /* length n*n row-major, owned by this struct */
  double *argument;  /* length n, owned by this struct */

  /* blather arrays, only allocated when opts.blather != 0.
   *
   * path_len is normally == iterations, but on a mid-loop objfun error
   * (error != 0 and error came from the objfun(theta.try) call rather
   * than the final post-loop re-evaluation) try_len == iterations - 1,
   * because the "pre" values (argpath/r/valpath) are recorded at the
   * *start* of each loop iteration while the "try" values are only
   * recorded once that iteration's trial point evaluates successfully.
   * This mirrors upstream trust()'s vector-growing behavior exactly. */
  int path_len;
  int try_len;
  double *argpath;    /* path_len x n, row-major */
  double *argtry;     /* try_len x n, row-major */
  int *steptype;      /* try_len; see trust_steptype_t below */
  int *accept;        /* try_len; 0/1 */
  double *r;          /* path_len */
  double *rho;        /* try_len */
  double *valpath;    /* path_len */
  double *valtry;     /* try_len */
  double *preddiff;   /* try_len */
  double *stepnorm;   /* try_len */
} trust_result_t;

typedef enum {
  TRUST_STEP_NEWTON = 0,
  TRUST_STEP_EASY_EASY = 1,
  TRUST_STEP_HARD_EASY = 2,
  TRUST_STEP_HARD_HARD = 3
} trust_steptype_t;

/* Zero-initializes a trust_result_t (does not allocate). */
void trust_result_zero(trust_result_t *res);

/* Frees every buffer owned by *res (safe to call on a zeroed or already
 * partially-filled result) and zeroes it. Uses free() only, matching the
 * malloc/calloc used internally -- safe to call from any thread. */
void trust_result_free(trust_result_t *res);

/* Thread-safe core solver: no R API calls, safe from OpenMP/parallel code. */
int trust_solve_c(int n, const double *parinit, trust_c_objfun_t objfun,
                   void *userdata, const trust_options_t *opts,
                   trust_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* TRUST_TYPES_H */
