// Native C objectives and OpenMP stress-test shims for the BOBYQA, NEWUOA
// and Steihaug ports, used only by tests/testthat (mirrors
// trust_test_objfuns.cpp's trust_openmp_stress_test for trust_solve_c()).
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <vector>
#include "minqa_types.h"
#include "steihaug_types.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {

// Extended Rosenbrock, sum_i 100 (x_{i+1} - x_i^2)^2 + (1 - x_i)^2.
double ext_rosen(int n, const double *x) {
  double f = 0.0;
  for (int i = 0; i + 1 < n; i++) {
    double t = x[i + 1] - x[i] * x[i];
    f += 100.0 * t * t + (1.0 - x[i]) * (1.0 - x[i]);
  }
  return f;
}

int minqa_test_rosen(int n, const double *x, double *f, void *) {
  *f = ext_rosen(n, x);
  return 0;
}

void ext_rosen_grad(int n, const double *x, double *g) {
  for (int i = 0; i < n; i++) g[i] = 0.0;
  for (int i = 0; i + 1 < n; i++) {
    double t = x[i + 1] - x[i] * x[i];
    g[i] += -400.0 * x[i] * t - 2.0 * (1.0 - x[i]);
    g[i + 1] += 200.0 * t;
  }
}

void ext_rosen_hess(int n, const double *x, double *h) {
  for (int i = 0; i < n * n; i++) h[i] = 0.0;
  for (int i = 0; i + 1 < n; i++) {
    h[i * n + i] += 1200.0 * x[i] * x[i] - 400.0 * x[i + 1] + 2.0;
    h[i * n + i + 1] += -400.0 * x[i];
    h[(i + 1) * n + i] += -400.0 * x[i];
    h[(i + 1) * n + i + 1] += 200.0;
  }
}

int steihaug_test_rosen(int n, const double *x, double *value, double *gradient,
                        double *hessian, void *) {
  if (value != nullptr) *value = ext_rosen(n, x);
  if (gradient != nullptr) ext_rosen_grad(n, x, gradient);
  if (hessian != nullptr) ext_rosen_hess(n, x, hessian);
  return 0;
}

int steihaug_test_rosen_hv(int n, const double *x, const double *v, double *hv,
                           void *) {
  std::vector<double> h(static_cast<size_t>(n) * n);
  ext_rosen_hess(n, x, h.data());
  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) s += h[i * n + j] * v[j];
    hv[i] = s;
  }
  return 0;
}

struct RunOut {
  std::vector<double> par;
  double value = 0.0;
  int evals = 0, code = 0;
};

// solver: 0 = bobyqa, 1 = newuoa, 2 = steihaug exact, 3 = steihaug matrix-free
RunOut run_one(int solver, const std::vector<double> &start) {
  int n = static_cast<int>(start.size());
  RunOut o;
  if (solver <= 1) {
    minqa_options_t opts = minqa_options_default(n, start.data());
    opts.npt = 2 * n + 1;
    minqa_result_t res;
    std::vector<double> lower(n, -10.0), upper(n, 10.0);
    if (solver == 0) {
      bobyqa_solve_c(n, start.data(), lower.data(), upper.data(),
                     minqa_test_rosen, nullptr, &opts, &res);
    } else {
      newuoa_solve_c(n, start.data(), minqa_test_rosen, nullptr, &opts, &res);
    }
    if (res.par != nullptr) o.par.assign(res.par, res.par + n);
    o.value = res.fval;
    o.evals = res.feval;
    o.code = res.ierr;
    minqa_result_free(&res);
  } else {
    steihaug_options_t opts = steihaug_options_default();
    steihaug_result_t res;
    steihaug_solve_c(n, start.data(), steihaug_test_rosen,
                     solver == 3 ? steihaug_test_rosen_hv : nullptr, nullptr,
                     &opts, &res);
    if (res.argument != nullptr) o.par.assign(res.argument, res.argument + n);
    o.value = res.value;
    o.evals = static_cast<int>(res.cost_evals + res.hessvec_evals);
    o.code = res.status;
    steihaug_result_free(&res);
  }
  return o;
}

}  // namespace

// Runs nStarts distinct extended-Rosenbrock problems once sequentially and
// once across OpenMP threads, and reports whether the results match
// bitwise. solver: "bobyqa", "newuoa", "steihaug", "steihaug_hessvec".
// [[Rcpp::export]]
List solver_openmp_stress_test(std::string solver, int n, int nStarts,
                               int nThreadsRequested) {
  int s = solver == "bobyqa" ? 0 : solver == "newuoa" ? 1
          : solver == "steihaug" ? 2 : solver == "steihaug_hessvec" ? 3 : -1;
  if (s < 0) Rcpp::stop("unknown solver");
  std::vector<std::vector<double> > starts(nStarts, std::vector<double>(n));
  for (int i = 0; i < nStarts; i++)
    for (int j = 0; j < n; j++)
      starts[i][j] = (j % 2 == 0 ? -1.2 : 1.0) + 0.013 * i - 0.004 * j;

  std::vector<RunOut> ref(nStarts), par(nStarts);
  for (int i = 0; i < nStarts; i++) ref[i] = run_one(s, starts[i]);

  std::vector<int> threadIds(nStarts, 0);
#ifdef _OPENMP
  int usedThreads = std::max(1, std::min(nThreadsRequested, omp_get_max_threads()));
#pragma omp parallel for num_threads(usedThreads) schedule(static)
#endif
  for (int i = 0; i < nStarts; i++) {
#ifdef _OPENMP
    threadIds[i] = omp_get_thread_num();
#endif
    par[i] = run_one(s, starts[i]);
  }

  bool identical = true;
  int nOk = 0;
  for (int i = 0; i < nStarts; i++) {
    if (ref[i].par != par[i].par || ref[i].evals != par[i].evals ||
        ref[i].code != par[i].code ||
        std::memcmp(&ref[i].value, &par[i].value, sizeof(double)) != 0)
      identical = false;
    bool near = !ref[i].par.empty();
    for (double v : ref[i].par) near = near && std::fabs(v - 1.0) < 1e-3;
    if (near) nOk++;
  }
  (void)nThreadsRequested;
  return List::create(
      _["identical"] = identical, _["nNearOptimum"] = nOk,
      _["distinctThreadsUsed"] =
          static_cast<int>(std::set<int>(threadIds.begin(), threadIds.end()).size()),
      _["openmpAvailable"] =
#ifdef _OPENMP
          true
#else
          false
#endif
  );
}

namespace {
int minqa_fail_after(int n, const double *x, double *f, void *ud) {
  int *left = static_cast<int *>(ud);
  if ((*left)-- <= 0) return -1;
  return minqa_test_rosen(n, x, f, nullptr);
}
int steihaug_fail_after(int n, const double *x, double *value, double *gradient,
                        double *hessian, void *ud) {
  int *left = static_cast<int *>(ud);
  if ((*left)-- <= 0) return -1;
  return steihaug_test_rosen(n, x, value, gradient, hessian, nullptr);
}
int steihaug_hv_fail(int, const double *, const double *, double *, void *) {
  return -1;
}
}  // namespace

// Return codes of the thread-safe C entry points on invalid input and
// failing callbacks (paths the R wrappers never reach).
// [[Rcpp::export]]
IntegerVector solver_c_api_edge_test() {
  double start[2] = {-1.2, 1.0}, zero[2] = {0.0, 0.0};
  double lower[2] = {0.0, 0.0}, upper[2] = {1e-8, 1.0};
  minqa_result_t mr;
  steihaug_result_t sr;
  IntegerVector out;
  auto push = [&](const char *nm, int code) {
    out.push_back(code, nm);
  };

  push("newuoa_null_opts", newuoa_solve_c(2, start, minqa_test_rosen, nullptr,
                                          nullptr, &mr));
  minqa_result_free(&mr);
  push("bobyqa_null_opts", bobyqa_solve_c(2, start, nullptr, nullptr,
                                          minqa_test_rosen, nullptr, nullptr,
                                          &mr));
  minqa_result_free(&mr);
  push("newuoa_zero_start", newuoa_solve_c(2, zero, minqa_test_rosen, nullptr,
                                           nullptr, &mr));
  push("bobyqa_negative_n", bobyqa_solve_c(-1, start, nullptr, nullptr,
                                           minqa_test_rosen, nullptr, nullptr,
                                           &mr));
  push("newuoa_null_result", newuoa_solve_c(2, start, minqa_test_rosen,
                                            nullptr, nullptr, nullptr));
  minqa_options_t mo = minqa_options_default(2, start);
  mo.npt = 2;
  push("newuoa_bad_npt", newuoa_solve_c(2, start, minqa_test_rosen, nullptr,
                                        &mo, &mr));
  minqa_result_free(&mr);
  push("bobyqa_bad_npt", bobyqa_solve_c(2, start, nullptr, nullptr,
                                        minqa_test_rosen, nullptr, &mo, &mr));
  minqa_result_free(&mr);
  mo = minqa_options_default(2, start);
  double tight[2] = {0.0, 0.5};
  push("bobyqa_range", bobyqa_solve_c(2, tight, lower, upper, minqa_test_rosen,
                                      nullptr, &mo, &mr));
  minqa_result_free(&mr);
  int left = 7;
  push("newuoa_objfun_error", newuoa_solve_c(2, start, minqa_fail_after, &left,
                                             &mo, &mr));
  minqa_result_free(&mr);
  left = 7;
  push("bobyqa_objfun_error", bobyqa_solve_c(2, start, nullptr, nullptr,
                                             minqa_fail_after, &left, &mo, &mr));
  minqa_result_free(&mr);
  // A non-finite start reaches calfun on the first evaluation, which aborts
  // (minqa stops with "non-finite x values not allowed in calfun").
  double nanstart[2] = {NAN, 1.0};
  int nanCode = newuoa_solve_c(2, nanstart, minqa_test_rosen, nullptr, &mo, &mr);
  minqa_result_free(&mr);
  push("newuoa_nonfinite_x", nanCode);

  steihaug_options_t so = steihaug_options_default();
  push("steihaug_ok", steihaug_solve_c(2, start, steihaug_test_rosen, nullptr,
                                       nullptr, &so, &sr));
  steihaug_result_free(&sr);
  push("steihaug_null_opts", steihaug_solve_c(2, start, steihaug_test_rosen,
                                              nullptr, nullptr, nullptr, &sr));
  steihaug_result_free(&sr);
  so.eta = 0.5;
  push("steihaug_bad_eta", steihaug_solve_c(2, start, steihaug_test_rosen,
                                            nullptr, nullptr, &so, &sr));
  steihaug_result_free(&sr);
  so = steihaug_options_default();
  so.abs_step_tol = NAN;
  push("steihaug_nan_tol", steihaug_solve_c(2, start, steihaug_test_rosen,
                                            nullptr, nullptr, &so, &sr));
  steihaug_result_free(&sr);
  so = steihaug_options_default();
  left = 5;
  push("steihaug_objfun_error", steihaug_solve_c(2, start, steihaug_fail_after,
                                                 nullptr, &left, &so, &sr));
  steihaug_result_free(&sr);
  push("steihaug_hessvec_error", steihaug_solve_c(2, start, steihaug_test_rosen,
                                                  steihaug_hv_fail, nullptr, &so,
                                                  &sr));
  steihaug_result_free(&sr);
  push("steihaug_zero_n", steihaug_solve_c(0, start, steihaug_test_rosen,
                                           nullptr, nullptr, &so, &sr));
  steihaug_result_free(&sr);
  return out;
}
