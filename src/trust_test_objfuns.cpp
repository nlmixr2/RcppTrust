// Native C objective functions and an Rcpp test shim for trust_solve_c(),
// used only while building the port incrementally (to validate the
// thread-safe C path directly against upstream trust::trust() output,
// and later for the OpenMP thread-safety stress test) -- not part of the
// public C interface itself.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include "trust_types.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Rosenbrock's function, matching tests/foo.R's objfun exactly:
// f = 100*(x2 - x1^2)^2 + (1 - x1)^2
extern "C" int trust_test_rosenbrock(int n, const double *par, double *value,
                                      double *gradient, double *hessian,
                                      void *userdata) {
  (void)userdata;
  if (n != 2) return -1;
  double x1 = par[0], x2 = par[1];
  double t = x2 - x1 * x1;
  *value = 100.0 * t * t + (1.0 - x1) * (1.0 - x1);
  gradient[0] = -400.0 * x1 * t - 2.0 * (1.0 - x1);
  gradient[1] = 200.0 * t;
  hessian[0] = -400.0 * x2 + 1200.0 * x1 * x1 + 2.0;  // d2f/dx1^2
  hessian[1] = -400.0 * x1;                            // d2f/dx1dx2
  hessian[2] = -400.0 * x1;                             // d2f/dx2dx1
  hessian[3] = 200.0;                                   // d2f/dx2^2
  return 0;
}

// [[Rcpp::export]]
List trust_solve_c_rosenbrock_test(NumericVector parinit, double rinit, double rmax,
                                    int iterlim, double fterm, double mterm,
                                    bool minimize, bool blather) {
  int n = parinit.size();
  trust_options_t opts = trust_options_default(rinit, rmax);
  opts.iterlim = iterlim;
  opts.fterm = fterm;
  opts.mterm = mterm;
  opts.minimize = minimize ? 1 : 0;
  opts.blather = blather ? 1 : 0;

  trust_result_t res;
  int rc = trust_solve_c(n, parinit.begin(), trust_test_rosenbrock, nullptr, &opts, &res);

  List out;
  out["rc"] = rc;
  out["error"] = res.error;
  out["converged"] = res.converged != 0;
  out["iterations"] = res.iterations;
  if (res.argument != nullptr)
    out["argument"] = NumericVector(res.argument, res.argument + n);
  if (res.gradient != nullptr) {
    out["value"] = res.value;
    out["gradient"] = NumericVector(res.gradient, res.gradient + n);
    NumericMatrix H(n, n);
    for (int i = 0; i < n * n; i++) H[i] = res.hessian[i];
    out["hessian"] = H;
  }
  if (blather && res.path_len > 0) {
    NumericMatrix argpath(res.path_len, n), argtry(res.try_len, n);
    for (int i = 0; i < res.path_len; i++)
      for (int j = 0; j < n; j++) argpath(i, j) = res.argpath[i * n + j];
    for (int i = 0; i < res.try_len; i++)
      for (int j = 0; j < n; j++) argtry(i, j) = res.argtry[i * n + j];
    out["argpath"] = argpath;
    out["argtry"] = argtry;
    out["r"] = NumericVector(res.r, res.r + res.path_len);
    out["valpath"] = NumericVector(res.valpath, res.valpath + res.path_len);
    out["rho"] = NumericVector(res.rho, res.rho + res.try_len);
    out["valtry"] = NumericVector(res.valtry, res.valtry + res.try_len);
    out["preddiff"] = NumericVector(res.preddiff, res.preddiff + res.try_len);
    out["stepnorm"] = NumericVector(res.stepnorm, res.stepnorm + res.try_len);
    CharacterVector steptype(res.try_len);
    static const char *names[] = {"Newton", "easy-easy", "hard-easy", "hard-hard"};
    for (int i = 0; i < res.try_len; i++) steptype[i] = names[res.steptype[i]];
    out["steptype"] = steptype;
    LogicalVector accept(res.try_len);
    for (int i = 0; i < res.try_len; i++) accept[i] = res.accept[i] != 0;
    out["accept"] = accept;
  }

  trust_result_free(&res);
  return out;
}

// Runs trust_solve_c() on nStarts distinct Rosenbrock starting points,
// once sequentially (the reference) and once with each start's solve
// dispatched to a different OpenMP thread, then compares the two sets of
// results. trust_solve_c() takes no locks and touches no shared mutable
// state (no R API calls, no globals -- each call's arma objects are
// thread-local automatic storage), so this is expected to reproduce the
// sequential results exactly regardless of thread count.
// [[Rcpp::export]]
List trust_openmp_stress_test(int nStarts, int nThreadsRequested) {
  std::vector<std::vector<double> > starts(nStarts, std::vector<double>(2));
  for (int i = 0; i < nStarts; i++) {
    starts[i][0] = 3.0 + 0.01 * i;
    starts[i][1] = 1.0 - 0.007 * i;
  }

  trust_options_t opts = trust_options_default(1.0, 5.0);

  std::vector<double> refArg(nStarts * 2), refVal(nStarts);
  std::vector<int> refConverged(nStarts), refIter(nStarts);
  for (int i = 0; i < nStarts; i++) {
    trust_result_t res;
    trust_solve_c(2, starts[i].data(), trust_test_rosenbrock, nullptr, &opts, &res);
    refArg[2 * i] = res.argument[0];
    refArg[2 * i + 1] = res.argument[1];
    refVal[i] = res.value;
    refConverged[i] = res.converged;
    refIter[i] = res.iterations;
    trust_result_free(&res);
  }

  std::vector<double> parArg(nStarts * 2), parVal(nStarts);
  std::vector<int> parConverged(nStarts), parIter(nStarts);
  std::vector<int> threadIds(nStarts, -1);

#ifdef _OPENMP
  int usedThreads = std::min(nThreadsRequested, omp_get_max_threads());
  if (usedThreads < 1) usedThreads = 1;
#pragma omp parallel for num_threads(usedThreads) schedule(static)
#endif
  for (int i = 0; i < nStarts; i++) {
#ifdef _OPENMP
    threadIds[i] = omp_get_thread_num();
#else
    threadIds[i] = 0;
#endif
    trust_result_t res;
    trust_solve_c(2, starts[i].data(), trust_test_rosenbrock, nullptr, &opts, &res);
    parArg[2 * i] = res.argument[0];
    parArg[2 * i + 1] = res.argument[1];
    parVal[i] = res.value;
    parConverged[i] = res.converged;
    parIter[i] = res.iterations;
    trust_result_free(&res);
  }

  double maxArgDiff = 0.0, maxValDiff = 0.0;
  bool allConvergedMatch = true, allIterMatch = true;
  for (int i = 0; i < nStarts; i++) {
    maxArgDiff = std::max(maxArgDiff, std::fabs(refArg[2 * i] - parArg[2 * i]));
    maxArgDiff = std::max(maxArgDiff, std::fabs(refArg[2 * i + 1] - parArg[2 * i + 1]));
    maxValDiff = std::max(maxValDiff, std::fabs(refVal[i] - parVal[i]));
    if (refConverged[i] != parConverged[i]) allConvergedMatch = false;
    if (refIter[i] != parIter[i]) allIterMatch = false;
  }

#ifdef _OPENMP
  int distinctThreadsUsed =
      static_cast<int>(std::set<int>(threadIds.begin(), threadIds.end()).size());
#else
  int distinctThreadsUsed = 1;
#endif

  return List::create(
      _["maxArgDiff"] = maxArgDiff, _["maxValDiff"] = maxValDiff,
      _["allConvergedMatch"] = allConvergedMatch, _["allIterMatch"] = allIterMatch,
      _["distinctThreadsUsed"] = distinctThreadsUsed,
      _["openmpAvailable"] =
#ifdef _OPENMP
          true
#else
          false
#endif
  );
}
