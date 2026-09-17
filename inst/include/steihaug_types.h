#ifndef STEIHAUG_TYPES_H
#define STEIHAUG_TYPES_H

/* Plain C types for the thread-safe port of the Steihaug (truncated
 * conjugate-gradient) trust-region Newton minimizer from the Rust crate
 * 'basin' (https://github.com/jolars/basin, MIT OR Apache-2.0, ported
 * from commit 86cc0d6: TrustRegion<Steihaug> / TrustRegion::matrix_free()
 * driven by basin's Executor). No R API and no C++ types appear here so
 * this header can be included from a consumer's plain C code as well. */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Objective function contract.
 *
 * n        : number of parameters
 * x        : parameter vector, length n (input)
 * value    : output objective value, or NULL when not requested
 * gradient : output buffer of length n, or NULL when not requested
 * hessian  : output buffer of length n*n, row-major (symmetric, so
 *            column-major storage is equivalent), or NULL when not
 *            requested
 * userdata : opaque pointer passed through unchanged
 *
 * The solver requests exactly what basin's problem calls request, one
 * callback per basin call, so evaluation counts match basin:
 *   - value + gradient      once at initialization (cost_and_gradient)
 *   - hessian only          once per outer iteration (exact mode only)
 *   - value only            once per trial step
 *   - gradient only         once after each accepted step
 * A non-finite value is allowed (it is treated as a rejected step, as in
 * basin). Return 0 on success, <0 to abort the run (basin propagates a
 * problem error as Err): the solver returns STEIHAUG_ERR_OBJFUN. */
typedef int (*steihaug_c_objfun_t)(int n, const double *x, double *value,
                                   double *gradient, double *hessian,
                                   void *userdata);

/* Optional Hessian-vector product hv = H(x) v (both length n). Passing a
 * non-NULL function selects basin's MatrixFree mode: the Hessian is never
 * requested from objfun. Return 0 on success, <0 to abort the run
 * (STEIHAUG_ERR_HESSVEC). */
typedef int (*steihaug_c_hessvec_t)(int n, const double *x, const double *v,
                                    double *hv, void *userdata);

/* Tolerance sentinel: a tolerance field set to any NEGATIVE value (use
 * STEIHAUG_DISABLED) is disabled. Enabled tolerances must be finite and
 * >= 0 (basin's optional_tolerance assertion); NaN or +Inf is rejected
 * with STEIHAUG_ERR_INVALID_OPTIONS. */
#define STEIHAUG_DISABLED (-1.0)

typedef struct {
  /* --- TrustRegion solver settings (basin defaults) --- */
  double initial_radius; /* Delta_0, default 1; must be > 0            */
  double max_radius;     /* Delta_max, default 100; must be > 0         */
  double eta;            /* acceptance threshold, default 0.125, [0,1/4) */
  int max_inner;         /* radius reductions per outer iteration,
                            default 10, >= 1                            */
  /* --- Steihaug CG settings --- */
  int cg_max_iter;       /* CG iteration cap, 0 = n (default), else >= 1 */
  double kappa;          /* forcing: stop when ||r|| < min(kappa,       */
  double theta;          /*   ||g||^theta) ||g||; defaults 0.5 / 0.5.
                            kappa finite in [0,1), theta finite >= 0     */

  /* --- Executor run controls (checked before every iteration, in this
   *     order: max_iter, max_cost_evals, max_gradient_evals, target) --- */
  int max_iter;             /* default 100 (basin's executor: 1000); >= 0 */
  long long max_cost_evals; /* < 0 disables (default)                    */
  long long max_gradient_evals; /* < 0 disables (default); counts only
                                   gradient evaluations, like basin      */
  int has_target_cost;      /* 0/1, default 0                            */
  double target_cost;       /* stop when best cost <= target; finite     */

  /* --- Solver convergence checks (after the run controls, in this order:
   *     absolute gradient, relative gradient, absolute step, relative
   *     step, absolute cost change, relative cost change). Negative =
   *     disabled. --- */
  double abs_gradient_tol;  /* ||g||_2 <= tol. Default 1e-8 (basin: off) */
  double rel_gradient_tol;  /* ||g||_2 <= tol ||g_0||_2, anchored on the
                               first finite gradient checked. Default off */
  double abs_step_tol;      /* ||x_k - x_{k-1}||_2 <= tol. Default off.
                               NB: an iteration whose inner attempts all
                               fail leaves x unchanged and fires this.  */
  double rel_step_tol;      /* ||x_k - x_{k-1}|| <= tol ||x_k||. Off      */
  double abs_cost_change_tol; /* |f_{k-1} - f_k| <= tol. Default off     */
  double rel_cost_change_tol; /* |f_{k-1} - f_k| <= tol |f_{k-1}|. Off   */

  int trace; /* 1 = record one trace row per subproblem attempt */
} steihaug_options_t;

static inline steihaug_options_t steihaug_options_default(void) {
  steihaug_options_t o;
  o.initial_radius = 1.0;
  o.max_radius = 100.0;
  o.eta = 0.125;
  o.max_inner = 10;
  o.cg_max_iter = 0;
  o.kappa = 0.5;
  o.theta = 0.5;
  o.max_iter = 100;
  o.max_cost_evals = -1;
  o.max_gradient_evals = -1;
  o.has_target_cost = 0;
  o.target_cost = 0.0;
  o.abs_gradient_tol = 1e-8;
  o.rel_gradient_tol = STEIHAUG_DISABLED;
  o.abs_step_tol = STEIHAUG_DISABLED;
  o.rel_step_tol = STEIHAUG_DISABLED;
  o.abs_cost_change_tol = STEIHAUG_DISABLED;
  o.rel_cost_change_tol = STEIHAUG_DISABLED;
  o.trace = 0;
  return o;
}

/* Termination reasons (>= 0, mirroring basin's TerminationReason variants
 * reachable by this solver) and errors (< 0). */
typedef enum {
  STEIHAUG_MAX_ITER = 0,
  STEIHAUG_MAX_COST_EVALS = 1,
  STEIHAUG_MAX_GRADIENT_EVALS = 2,
  STEIHAUG_TARGET_COST = 3,
  STEIHAUG_GRADIENT_TOLERANCE = 4,
  STEIHAUG_RELATIVE_GRADIENT_TOLERANCE = 5,
  STEIHAUG_PARAM_TOLERANCE = 6,          /* absolute step tolerance */
  STEIHAUG_RELATIVE_PARAM_TOLERANCE = 7, /* relative step tolerance */
  STEIHAUG_COST_TOLERANCE = 8,           /* absolute cost change    */
  STEIHAUG_RELATIVE_COST_TOLERANCE = 9,  /* relative cost change    */
  STEIHAUG_SOLVER_CONVERGED = 10, /* predicted reduction <= 0 (g == 0)  */
  STEIHAUG_SOLVER_FAILED = 11,    /* non-finite predicted reduction      */

  STEIHAUG_ERR_OBJFUN = -1,          /* objfun returned < 0              */
  STEIHAUG_ERR_HESSVEC = -2,         /* hessvec returned < 0             */
  STEIHAUG_ERR_INVALID_OPTIONS = -3, /* n < 1, NULL pointers, bad opts   */
  STEIHAUG_ERR_NOMEM = -4,           /* allocation failure               */
  STEIHAUG_ERR_INTERNAL = -99
} steihaug_status_t;

typedef struct {
  int n;
  int status;     /* steihaug_status_t */
  int converged;  /* 1 for the tolerance / target / SolverConverged
                     reasons (4..10 and 3), else 0 */
  int iterations; /* basin's state.iter: completed outer iterations
                     (a SolverConverged/SolverFailed stop does not count) */

  /* Final state (basin's state.param / state.cost / state.gradient). On
   * an ERR_OBJFUN/ERR_HESSVEC abort these hold the state at the start of
   * the failing outer iteration (basin itself discards the state). */
  double value;
  double *argument; /* length n, malloc-owned */
  double *gradient; /* length n, malloc-owned */
  double radius;    /* trust radius after the run */

  /* Evaluation counts, incremented before each call like basin's Problem
   * wrapper (so a failing call is counted). */
  long long cost_evals;     /* calls that requested value */
  long long gradient_evals; /* calls that requested gradient */
  long long hessian_evals;
  long long hessvec_evals;

  /* Trace (trace = 1): one row per subproblem attempt, length trace_len.
   * All arrays malloc-owned. */
  int trace_len;
  int *trace_iter;          /* outer iteration index (state.iter)         */
  double *trace_radius;     /* radius used for this attempt               */
  double *trace_cost;       /* cost at the current iterate                */
  double *trace_gradnorm;   /* ||g||_2 at the current iterate             */
  double *trace_trial_cost; /* cost at x + d (NaN if not evaluated)       */
  double *trace_pred;       /* predicted reduction m(0) - m(d)            */
  double *trace_rho;        /* reduction ratio (NaN if not evaluated)     */
  double *trace_step_norm;  /* ||d||_2 (NaN if not evaluated)             */
  int *trace_accept;        /* 0/1                                        */
  int *trace_hit_boundary;  /* 0/1                                        */
  int *trace_cg_iter;       /* CG iterations performed in this attempt    */
} steihaug_result_t;

/* Zero-initializes a result (does not allocate). */
void steihaug_result_zero(steihaug_result_t *res);

/* Frees every malloc-owned buffer and zeroes the struct. Safe on a zeroed
 * or partially filled result, from any thread. */
void steihaug_result_free(steihaug_result_t *res);

/* Thread-safe solver: no R API, no globals or statics; all buffers are
 * owned by the call. `result` is zeroed on entry (free any previous
 * contents first). Returns result->status. */
int steihaug_solve_c(int n, const double *parinit, steihaug_c_objfun_t objfun,
                     steihaug_c_hessvec_t hessvec_or_null, void *userdata,
                     const steihaug_options_t *opts,
                     steihaug_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* STEIHAUG_TYPES_H */
