// Thread-safe C++17 port of basin's Steihaug truncated-CG trust-region
// Newton minimizer (https://github.com/jolars/basin, MIT OR Apache-2.0,
// commit 86cc0d6, crate `basin` v1.12.0).
//
// Ported sources (paths relative to crates/basin/src):
//   solver/trust_region.rs           tr_init, tr_next_iter, both modes,
//                                    tau_to_boundary, model_decrease_from_bd
//   solver/trust_region/steihaug.rs  Steihaug::solve_with
//   core/executor.rs                 Executor::into_stepper, step_once
//   core/run_control.rs              RunControl::check
//   core/convergence.rs              ConfiguredSolver::check_convergence,
//                                    GradientChecks, StepChecks, CostChecks
//   core/termination.rs              GradientTolerance,
//                                    RelativeGradientTolerance,
//                                    CostTolerance, RelativeCostTolerance,
//                                    TargetCost
//   core/math/vec.rs, dense.rs       Vec<f64> / DenseMatrix kernels
//
// The floating-point kernels reproduce the Vec<f64>/DenseMatrix backend's
// order of operations (left-to-right sums seeded with -0.0, as Rust's
// `Iterator::sum` for f64 is) so iterates can match basin bit-for-bit.
// FMA contraction is disabled below via pragmas (rustc never contracts).
//
// Thread safety: no statics, no globals, no R API. Every buffer is owned
// by the call; callbacks receive the caller's userdata.

#include "steihaug_types.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <new>
#include <vector>

// Forbid FMA contraction (a*b + c -> fma) in this translation unit so the
// kernels below keep basin's (rustc's) exact IEEE operation sequence even
// when the package is built with -march=native / -mfma. Placed after the
// system headers so only this file's own functions are affected.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

namespace {

using Vec = std::vector<double>;

// ---------------------------------------------------------------------------
// core/math/vec.rs and core/math/dense.rs kernels
// ---------------------------------------------------------------------------

// vec.rs `impl Dot for Vec<F>`: iter().zip().map(a*b).sum() (seed -0.0).
inline double dot(const Vec &a, const Vec &b) {
  double s = -0.0;
  const size_t n = a.size();
  for (size_t i = 0; i < n; i++) s = s + a[i] * b[i];
  return s;
}

// vec.rs `impl NormSquared for Vec<F>`: iter().map(x*x).sum().
inline double norm_squared(const Vec &a) {
  double s = -0.0;
  const size_t n = a.size();
  for (size_t i = 0; i < n; i++) s = s + a[i] * a[i];
  return s;
}

// vec.rs `impl ScaledAdd for Vec<F>`: x = x + scalar * y.
inline void scaled_add(Vec &x, double scalar, const Vec &y) {
  const size_t n = x.size();
  for (size_t i = 0; i < n; i++) x[i] = x[i] + scalar * y[i];
}

// vec.rs `impl ScaleInPlace for Vec<F>`: x = x * scalar.
inline void scale_in_place(Vec &x, double scalar) {
  for (double &v : x) v = v * scalar;
}

// vec.rs `impl NegInPlace for Vec<F>`: x = -x.
inline void neg_in_place(Vec &x) {
  for (double &v : x) v = -v;
}

// dense.rs `impl MatVec<Vec<F>> for DenseMatrix<F>`: row-major
// y[i] = sum_j A[i*cols + j] * x[j] (seed -0.0).
inline void matvec(const Vec &a, int n, const Vec &x, Vec &y) {
  const size_t nn = static_cast<size_t>(n);
  for (size_t i = 0; i < nn; i++) {
    const double *row = a.data() + i * nn;
    double s = -0.0;
    for (size_t j = 0; j < nn; j++) s = s + row[j] * x[j];
    y[i] = s;
  }
}

inline bool is_finite(double x) { return std::isfinite(x); }

// ---------------------------------------------------------------------------
// Problem wrapper: core/problem.rs `Problem<P>` counted calls
// ---------------------------------------------------------------------------

enum class Err { NONE = 0, OBJFUN, HESSVEC };

struct Problem {
  int n;
  steihaug_c_objfun_t objfun;
  steihaug_c_hessvec_t hessvec;
  void *userdata;
  long long cost_evals = 0;
  long long gradient_evals = 0;
  long long hessian_evals = 0;
  long long hessvec_evals = 0;

  // problem.rs Problem::cost_and_gradient (counts bumped before the call)
  bool cost_and_gradient(const Vec &x, double &cost, Vec &grad) {
    cost_evals += 1;
    gradient_evals += 1;
    grad.assign(static_cast<size_t>(n), 0.0);
    return objfun(n, x.data(), &cost, grad.data(), nullptr, userdata) >= 0;
  }
  // problem.rs Problem::cost
  bool cost(const Vec &x, double &cost) {
    cost_evals += 1;
    return objfun(n, x.data(), &cost, nullptr, nullptr, userdata) >= 0;
  }
  // problem.rs Problem::gradient
  bool gradient(const Vec &x, Vec &grad) {
    gradient_evals += 1;
    grad.assign(static_cast<size_t>(n), 0.0);
    return objfun(n, x.data(), nullptr, grad.data(), nullptr, userdata) >= 0;
  }
  // problem.rs Problem::hessian
  bool hessian(const Vec &x, Vec &h) {
    hessian_evals += 1;
    h.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    return objfun(n, x.data(), nullptr, nullptr, h.data(), userdata) >= 0;
  }
  // problem.rs Problem::hessian_product
  bool hessian_product(const Vec &x, const Vec &v, Vec &hv) {
    hessvec_evals += 1;
    hv.assign(static_cast<size_t>(n), 0.0);
    return hessvec(n, x.data(), v.data(), hv.data(), userdata) >= 0;
  }
};

// ---------------------------------------------------------------------------
// solver/trust_region.rs
// ---------------------------------------------------------------------------

// trust_region.rs `struct Step`.
struct Step {
  Vec d;
  double predicted_reduction;
  bool hit_boundary;
  int cg_iter;  // port-only diagnostic for the trace
};

// trust_region.rs `model_decrease_from_bd`: -g.dot(d) - half * d.dot(bd).
inline double model_decrease_from_bd(const Vec &g, const Vec &d,
                                     const Vec &bd) {
  const double half = 0.5;
  return -dot(g, d) - half * dot(d, bd);
}

// trust_region.rs `tau_to_boundary`.
inline double tau_to_boundary(const Vec &z, const Vec &d, double radius) {
  double dd = dot(d, d);
  double zd = dot(z, d);
  double zz = dot(z, z);
  double rr = radius * radius;
  double disc = zd * zd - dd * (zz - rr);
  disc = disc < 0.0 ? 0.0 : disc;
  return (-zd + std::sqrt(disc)) / dd;
}

// The Hessian-vector product closure `bv` of both modes: exact mode wraps
// DenseMatrix::matvec (trust_region.rs Subproblem for Steihaug, infallible),
// matrix-free mode wraps the counted Problem::hessian_product.
struct HessOp {
  Problem *problem;
  const Vec *x;
  const Vec *hessian;  // non-null in exact mode
  // Returns false on a hessvec failure (Err propagation).
  bool operator()(const Vec &v, Vec &out) {
    if (hessian != nullptr) {
      out.resize(v.size());
      matvec(*hessian, problem->n, v, out);
      return true;
    }
    return problem->hessian_product(*x, v, out);
  }
};

// steihaug.rs `Steihaug::solve_with`. Returns false on a product error.
bool steihaug_solve_with(const Vec &g, double radius, int cfg_max_iter,
                         double kappa, double theta, HessOp &bv, Step &out) {
  const int n = static_cast<int>(g.size());
  const int max_iter = cfg_max_iter > 0 ? cfg_max_iter : (n > 1 ? n : 1);
  out.cg_iter = 0;

  // z0 = 0, r0 = g, d0 = -r0.
  Vec z = g;
  scale_in_place(z, 0.0);
  Vec r = g;
  double r_dot = dot(r, r);

  const double g_norm = std::sqrt(r_dot);
  if (g_norm == 0.0) {
    out.d = std::move(z);
    out.predicted_reduction = 0.0;
    out.hit_boundary = false;
    return true;
  }

  const double half = 0.5;
  double power;
  if (theta == 0.0) {
    power = 1.0;
  } else if (theta == half) {
    power = std::sqrt(g_norm);
  } else if (theta == 1.0) {
    power = g_norm;
  } else {
    power = std::pow(g_norm, theta);
  }
  const double tol = (power < kappa ? power : kappa) * g_norm;

  Vec d = r;
  neg_in_place(d);
  Vec bd, bz, z_next, d_next;

  for (int it = 0; it < max_iter; it++) {
    out.cg_iter = it + 1;
    if (!bv(d, bd)) return false;
    double dbd = dot(d, bd);

    // Non-positive curvature: walk to the boundary and stop.
    if (dbd <= 0.0) {
      double tau = tau_to_boundary(z, d, radius);
      scaled_add(z, tau, d);
      if (!bv(z, bz)) return false;
      out.predicted_reduction = model_decrease_from_bd(g, z, bz);
      out.d = std::move(z);
      out.hit_boundary = true;
      return true;
    }

    double alpha = r_dot / dbd;

    // Tentative iterate z + alpha d; clip to the boundary if outside.
    z_next = z;
    scaled_add(z_next, alpha, d);
    if (std::sqrt(norm_squared(z_next)) >= radius) {
      double tau = tau_to_boundary(z, d, radius);
      scaled_add(z, tau, d);
      if (!bv(z, bz)) return false;
      out.predicted_reduction = model_decrease_from_bd(g, z, bz);
      out.d = std::move(z);
      out.hit_boundary = true;
      return true;
    }
    z.swap(z_next);

    // r <- r + alpha B d.
    scaled_add(r, alpha, bd);
    double r_dot_next = dot(r, r);
    double residual_norm = std::sqrt(r_dot_next);
    bool converged = residual_norm == 0.0 || residual_norm < tol;
    if (converged) {
      if (!bv(z, bz)) return false;
      out.predicted_reduction = model_decrease_from_bd(g, z, bz);
      out.d = std::move(z);
      out.hit_boundary = false;
      return true;
    }

    // d <- -r + beta d.
    double beta = r_dot_next / r_dot;
    d_next = r;
    neg_in_place(d_next);
    scaled_add(d_next, beta, d);
    d.swap(d_next);
    r_dot = r_dot_next;
  }

  // Iteration cap reached: best interior iterate.
  if (!bv(z, bz)) return false;
  out.predicted_reduction = model_decrease_from_bd(g, z, bz);
  out.d = std::move(z);
  out.hit_boundary = false;
  return true;
}

// ---------------------------------------------------------------------------
// State: core/state.rs BasicState<Vec<f64>, f64> (fields used here)
// ---------------------------------------------------------------------------

struct State {
  Vec param;
  double cost = 0.0;
  Vec gradient;
  long long iter = 0;
  bool has_best = false;  // best_param.is_some()
  double best_cost = HUGE_VAL;

  // state.rs BasicState::update_best (cost is always Some after init).
  void update_best() {
    if (!has_best || cost < best_cost) {
      has_best = true;
      best_cost = cost;
    }
  }
};

// ---------------------------------------------------------------------------
// Trace buffer (port-only diagnostic)
// ---------------------------------------------------------------------------

struct TraceRow {
  int iter;
  double radius, cost, gradnorm, trial_cost, pred, rho, step_norm;
  int accept, hit_boundary, cg_iter;
};

// ---------------------------------------------------------------------------
// convergence.rs: ConfiguredSolver checks (G = GradientChecks,
// X = StepChecks, C = CostChecks; simplex slot unused)
// ---------------------------------------------------------------------------

struct Convergence {
  // GradientChecks
  bool abs_grad_on, rel_grad_on;
  double abs_grad, rel_grad;
  bool rel_grad_anchor_set = false;
  double rel_grad_initial = 0.0;
  // StepChecks
  bool abs_step_on, rel_step_on;
  double abs_step, rel_step;
  bool has_last_param = false;
  Vec last_param;
  // CostChecks
  bool abs_cost_on, rel_cost_on;
  double abs_cost, rel_cost;
  bool abs_cost_has_last = false, rel_cost_has_last = false;
  double abs_cost_last = 0.0, rel_cost_last = 0.0;

  // convergence.rs `impl Check for GradientChecks`, with
  // termination.rs GradientTolerance / RelativeGradientTolerance.
  int check_gradient(const State &s) {
    if (!abs_grad_on && !rel_grad_on) return -1;  // `()` slot
    double ns = norm_squared(s.gradient);
    if (!is_finite(ns)) return -1;
    if (abs_grad_on) {
      if (norm_squared(s.gradient) <= abs_grad * abs_grad)
        return STEIHAUG_GRADIENT_TOLERANCE;
    }
    if (rel_grad_on) {
      double norm_sq = norm_squared(s.gradient);
      if (!rel_grad_anchor_set) {
        rel_grad_anchor_set = true;
        rel_grad_initial = norm_sq;
      }
      if (norm_sq <= rel_grad * rel_grad * rel_grad_initial)
        return STEIHAUG_RELATIVE_GRADIENT_TOLERANCE;
    }
    return -1;
  }

  // convergence.rs `impl Check for StepChecks`.
  int check_step(const State &s) {
    if (!abs_step_on && !rel_step_on) return -1;
    const Vec &current = s.param;
    // let last = self.last.replace(current.clone())?;
    if (!has_last_param) {
      has_last_param = true;
      last_param = current;
      return -1;
    }
    Vec last = last_param;
    last_param = current;
    Vec difference = current;
    scaled_add(difference, -1.0, last);
    double step = std::sqrt(norm_squared(difference));
    if (!is_finite(step)) return -1;
    if (abs_step_on && step <= abs_step) return STEIHAUG_PARAM_TOLERANCE;
    if (rel_step_on) {
      bool hit;
      if (rel_step == 0.0) {
        hit = step == 0.0;
      } else {
        double norm = std::sqrt(norm_squared(current));
        double bound;
        if (is_finite(norm)) {
          bound = rel_step * norm;
        } else {
          Vec scaled = current;
          scaled_add(scaled, -1.0, current);
          scaled_add(scaled, rel_step, current);
          bound = std::sqrt(norm_squared(scaled));
        }
        hit = step <= bound;
      }
      if (hit) return STEIHAUG_RELATIVE_PARAM_TOLERANCE;
    }
    return -1;
  }

  // convergence.rs `impl Check for CostChecks`, with termination.rs
  // CostTolerance / RelativeCostTolerance (both evaluated, absolute wins).
  int check_cost(const State &s) {
    if (!abs_cost_on && !rel_cost_on) return -1;
    const double curr = s.cost;
    if (!is_finite(curr)) {
      abs_cost_has_last = false;
      rel_cost_has_last = false;
      return -1;
    }
    int absolute = -1, relative = -1;
    if (abs_cost_on) {
      bool trig = abs_cost_has_last &&
                  (std::fabs(abs_cost_last - curr) <= abs_cost &&
                   is_finite(curr));
      abs_cost_has_last = true;
      abs_cost_last = curr;
      if (trig) absolute = STEIHAUG_COST_TOLERANCE;
    }
    if (rel_cost_on) {
      bool trig = rel_cost_has_last &&
                  (is_finite(curr) &&
                   std::fabs(rel_cost_last - curr) <=
                       rel_cost * std::fabs(rel_cost_last));
      rel_cost_has_last = true;
      rel_cost_last = curr;
      if (trig) relative = STEIHAUG_RELATIVE_COST_TOLERANCE;
    }
    return absolute >= 0 ? absolute : relative;
  }

  // convergence.rs ConfiguredSolver::check_convergence: gradient, then
  // step, then cost (the per-boundary memo never hits within one run:
  // every executor step either increments `iter` or stops).
  int check(const State &s) {
    int r = check_gradient(s);
    if (r >= 0) return r;
    r = check_step(s);
    if (r >= 0) return r;
    return check_cost(s);
  }
};

// Outcome of one `Solver::next_iter`.
enum class IterOutcome { CONTINUE, CONVERGED, FAILED, ERROR };

struct Solver {
  const steihaug_options_t *opts;
  int n;
  double radius;
  bool exact;
  Err err = Err::NONE;
  bool trace;
  std::vector<TraceRow> rows;

  // trust_region.rs tr_init.
  bool init(Problem &problem, State &state) {
    radius = opts->initial_radius;
    double cost = 0.0;
    Vec grad;
    if (!problem.cost_and_gradient(state.param, cost, grad)) {
      err = Err::OBJFUN;
      return false;
    }
    state.cost = cost;
    state.gradient = std::move(grad);
    return true;
  }

  // trust_region.rs tr_next_iter plus the mode-specific `next_iter`
  // bodies (`impl Solver for TrustRegion<Sub, F, ExactHessian>` evaluates
  // one Hessian up front; `MatrixFree` goes through hessian_product).
  // On success `next` holds the new state; on ERROR it is unspecified.
  IterOutcome next_iter(Problem &problem, const State &state, State &next) {
    Vec b;
    if (exact) {
      if (!problem.hessian(state.param, b)) {
        err = Err::OBJFUN;
        return IterOutcome::ERROR;
      }
    }
    next = state;
    const Vec &g = state.gradient;
    const double cost_old = state.cost;
    const double quarter = 0.25;
    const double three_quarters = 0.75;
    const double two = 2.0;

    HessOp bv{&problem, &state.param, exact ? &b : nullptr};
    Step step;
    Vec trial;
    const double gradnorm = trace ? std::sqrt(norm_squared(g)) : 0.0;

    for (int k = 0; k < opts->max_inner; k++) {
      const double radius_used = radius;
      if (!steihaug_solve_with(g, radius, opts->cg_max_iter, opts->kappa,
                               opts->theta, bv, step)) {
        err = Err::HESSVEC;
        return IterOutcome::ERROR;
      }
      TraceRow row{static_cast<int>(state.iter), radius_used, cost_old,
                   gradnorm, NAN, step.predicted_reduction, NAN, NAN, 0,
                   step.hit_boundary ? 1 : 0, step.cg_iter};

      if (!is_finite(step.predicted_reduction)) {
        if (trace) rows.push_back(row);
        return IterOutcome::FAILED;
      }
      if (step.predicted_reduction <= 0.0) {
        if (trace) rows.push_back(row);
        return IterOutcome::CONVERGED;
      }

      trial = state.param;
      scaled_add(trial, 1.0, step.d);
      double cost_trial = 0.0;
      if (!problem.cost(trial, cost_trial)) {
        err = Err::OBJFUN;
        return IterOutcome::ERROR;
      }

      double rho = (cost_old - cost_trial) / step.predicted_reduction;
      double step_norm = std::sqrt(norm_squared(step.d));

      if (rho < quarter || !is_finite(rho)) {
        radius = quarter * step_norm;
      } else if (rho > three_quarters && step.hit_boundary) {
        double grown = two * radius;
        radius = grown < opts->max_radius ? grown : opts->max_radius;
      }

      row.trial_cost = cost_trial;
      row.rho = rho;
      row.step_norm = step_norm;

      if (rho > opts->eta) {
        row.accept = 1;
        if (trace) rows.push_back(row);
        next.param = trial;
        next.cost = cost_trial;
        if (!problem.gradient(next.param, next.gradient)) {
          err = Err::OBJFUN;
          return IterOutcome::ERROR;
        }
        return IterOutcome::CONTINUE;
      }
      if (trace) rows.push_back(row);
    }
    // Inner attempts exhausted: keep iterate, shrunken radius.
    return IterOutcome::CONTINUE;
  }
};

bool valid_tol(double t) {
  // Negative = disabled; otherwise basin's optional_tolerance assertion.
  if (t < 0.0) return true;
  return is_finite(t) && t >= 0.0;
}

bool validate(int n, const double *parinit, steihaug_c_objfun_t objfun,
              const steihaug_options_t *o) {
  if (n < 1 || parinit == nullptr || objfun == nullptr || o == nullptr)
    return false;
  // TrustRegion::with_radius / with_max_radius: `> 0` (NaN rejected).
  if (!(o->initial_radius > 0.0)) return false;
  if (!(o->max_radius > 0.0)) return false;
  // with_eta: eta in [0, 1/4).
  if (!(o->eta >= 0.0 && o->eta < 0.25)) return false;
  // with_max_inner_attempts: >= 1.
  if (o->max_inner < 1) return false;
  // Steihaug::with_max_iter: >= 1 (0 here selects the default n).
  if (o->cg_max_iter < 0) return false;
  // with_forcing_parameters.
  if (!(is_finite(o->kappa) && o->kappa >= 0.0 && o->kappa < 1.0))
    return false;
  if (!(is_finite(o->theta) && o->theta >= 0.0)) return false;
  if (o->max_iter < 0) return false;
  // RunControl::target_cost: finite.
  if (o->has_target_cost && !is_finite(o->target_cost)) return false;
  if (!valid_tol(o->abs_gradient_tol) || !valid_tol(o->rel_gradient_tol) ||
      !valid_tol(o->abs_step_tol) || !valid_tol(o->rel_step_tol) ||
      !valid_tol(o->abs_cost_change_tol) || !valid_tol(o->rel_cost_change_tol))
    return false;
  // NaN tolerances: `t < 0` is false and is_finite fails -> rejected above.
  return true;
}

double *dup_doubles(const double *src, size_t count) {
  size_t bytes = (count > 0 ? count : 1) * sizeof(double);
  double *p = static_cast<double *>(malloc(bytes));
  if (p == nullptr) throw std::bad_alloc();
  if (src != nullptr && count > 0) memcpy(p, src, count * sizeof(double));
  return p;
}

int *alloc_ints(size_t count) {
  int *p = static_cast<int *>(malloc((count > 0 ? count : 1) * sizeof(int)));
  if (p == nullptr) throw std::bad_alloc();
  return p;
}

void fill_state(steihaug_result_t *res, const State &s, int n) {
  res->value = s.cost;
  res->iterations = static_cast<int>(s.iter);
  res->argument = dup_doubles(s.param.data(), static_cast<size_t>(n));
  res->gradient = dup_doubles(
      s.gradient.size() == static_cast<size_t>(n) ? s.gradient.data() : nullptr,
      static_cast<size_t>(n));
  if (s.gradient.size() != static_cast<size_t>(n))
    for (int i = 0; i < n; i++) res->gradient[i] = NAN;
}

void fill_trace(steihaug_result_t *res, const std::vector<TraceRow> &rows) {
  size_t m = rows.size();
  res->trace_iter = alloc_ints(m);
  res->trace_radius = dup_doubles(nullptr, m);
  res->trace_cost = dup_doubles(nullptr, m);
  res->trace_gradnorm = dup_doubles(nullptr, m);
  res->trace_trial_cost = dup_doubles(nullptr, m);
  res->trace_pred = dup_doubles(nullptr, m);
  res->trace_rho = dup_doubles(nullptr, m);
  res->trace_step_norm = dup_doubles(nullptr, m);
  res->trace_accept = alloc_ints(m);
  res->trace_hit_boundary = alloc_ints(m);
  res->trace_cg_iter = alloc_ints(m);
  for (size_t i = 0; i < m; i++) {
    const TraceRow &r = rows[i];
    res->trace_iter[i] = r.iter;
    res->trace_radius[i] = r.radius;
    res->trace_cost[i] = r.cost;
    res->trace_gradnorm[i] = r.gradnorm;
    res->trace_trial_cost[i] = r.trial_cost;
    res->trace_pred[i] = r.pred;
    res->trace_rho[i] = r.rho;
    res->trace_step_norm[i] = r.step_norm;
    res->trace_accept[i] = r.accept;
    res->trace_hit_boundary[i] = r.hit_boundary;
    res->trace_cg_iter[i] = r.cg_iter;
  }
  res->trace_len = static_cast<int>(m);
}

// executor.rs Executor::into_stepper + Stepper::run_to_end + step_once,
// with run_control.rs RunControl::check.
int solve(int n, const double *parinit, steihaug_c_objfun_t objfun,
          steihaug_c_hessvec_t hessvec, void *userdata,
          const steihaug_options_t *o, steihaug_result_t *res) {
  Problem problem{n, objfun, hessvec, userdata, 0, 0, 0, 0};
  Solver solver{o, n, o->initial_radius, hessvec == nullptr, Err::NONE,
                o->trace != 0, {}};

  Convergence conv;
  conv.abs_grad_on = o->abs_gradient_tol >= 0.0;
  conv.abs_grad = o->abs_gradient_tol;
  conv.rel_grad_on = o->rel_gradient_tol >= 0.0;
  conv.rel_grad = o->rel_gradient_tol;
  conv.abs_step_on = o->abs_step_tol >= 0.0;
  conv.abs_step = o->abs_step_tol;
  conv.rel_step_on = o->rel_step_tol >= 0.0;
  conv.rel_step = o->rel_step_tol;
  conv.abs_cost_on = o->abs_cost_change_tol >= 0.0;
  conv.abs_cost = o->abs_cost_change_tol;
  conv.rel_cost_on = o->rel_cost_change_tol >= 0.0;
  conv.rel_cost = o->rel_cost_change_tol;

  State state;
  state.param.assign(parinit, parinit + n);

  int status = STEIHAUG_ERR_INTERNAL;
  State next;
  bool ok = solver.init(problem, state);
  if (!ok) {
    // basin returns Err from into_stepper; report the starting point.
    status = STEIHAUG_ERR_OBJFUN;
    state.cost = NAN;
  } else {
    state.update_best();
    for (;;) {
      // RunControl::check.
      if (state.iter >= o->max_iter) {
        status = STEIHAUG_MAX_ITER;
        break;
      }
      if (o->max_cost_evals >= 0 && problem.cost_evals >= o->max_cost_evals) {
        status = STEIHAUG_MAX_COST_EVALS;
        break;
      }
      if (o->max_gradient_evals >= 0 &&
          problem.gradient_evals >= o->max_gradient_evals) {
        status = STEIHAUG_MAX_GRADIENT_EVALS;
        break;
      }
      if (o->has_target_cost && state.best_cost <= o->target_cost) {
        status = STEIHAUG_TARGET_COST;
        break;
      }
      // Solver::check_convergence.
      int reason = conv.check(state);
      if (reason >= 0) {
        status = reason;
        break;
      }
      // Solver::next_iter.
      IterOutcome out = solver.next_iter(problem, state, next);
      if (out == IterOutcome::ERROR) {
        status = solver.err == Err::HESSVEC ? STEIHAUG_ERR_HESSVEC
                                             : STEIHAUG_ERR_OBJFUN;
        break;
      }
      if (out == IterOutcome::CONVERGED || out == IterOutcome::FAILED) {
        // mid-iteration stop: state unchanged apart from update_best.
        state.update_best();
        status = out == IterOutcome::CONVERGED ? STEIHAUG_SOLVER_CONVERGED
                                               : STEIHAUG_SOLVER_FAILED;
        break;
      }
      state.param.swap(next.param);
      state.cost = next.cost;
      state.gradient.swap(next.gradient);
      state.iter += 1;
      state.update_best();
    }
  }

  res->n = n;
  res->status = status;
  res->converged = (status >= STEIHAUG_TARGET_COST &&
                    status <= STEIHAUG_SOLVER_CONVERGED)
                       ? 1
                       : 0;
  res->radius = solver.radius;
  res->cost_evals = problem.cost_evals;
  res->gradient_evals = problem.gradient_evals;
  res->hessian_evals = problem.hessian_evals;
  res->hessvec_evals = problem.hessvec_evals;
  fill_state(res, state, n);
  if (solver.trace) fill_trace(res, solver.rows);
  return status;
}

}  // namespace

extern "C" void steihaug_result_zero(steihaug_result_t *res) {
  if (res == nullptr) return;
  memset(res, 0, sizeof(*res));
}

extern "C" void steihaug_result_free(steihaug_result_t *res) {
  if (res == nullptr) return;
  free(res->argument);
  free(res->gradient);
  free(res->trace_iter);
  free(res->trace_radius);
  free(res->trace_cost);
  free(res->trace_gradnorm);
  free(res->trace_trial_cost);
  free(res->trace_pred);
  free(res->trace_rho);
  free(res->trace_step_norm);
  free(res->trace_accept);
  free(res->trace_hit_boundary);
  free(res->trace_cg_iter);
  memset(res, 0, sizeof(*res));
}

extern "C" int steihaug_solve_c(int n, const double *parinit,
                                steihaug_c_objfun_t objfun,
                                steihaug_c_hessvec_t hessvec_or_null,
                                void *userdata, const steihaug_options_t *opts,
                                steihaug_result_t *result) {
  if (result == nullptr) return STEIHAUG_ERR_INVALID_OPTIONS;
  steihaug_result_zero(result);
  if (!validate(n, parinit, objfun, opts)) {
    result->n = n;
    result->status = STEIHAUG_ERR_INVALID_OPTIONS;
    return result->status;
  }
  try {
    return solve(n, parinit, objfun, hessvec_or_null, userdata, opts, result);
  } catch (const std::bad_alloc &) {
    steihaug_result_free(result);
    result->n = n;
    result->status = STEIHAUG_ERR_NOMEM;
  } catch (...) {
    steihaug_result_free(result);
    result->n = n;
    result->status = STEIHAUG_ERR_INTERNAL;
  }
  return result->status;
}
