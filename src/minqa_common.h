#ifndef MINQA_COMMON_H
#define MINQA_COMMON_H

/* Shared, thread-safe replacement for minqa's R-bound glue (calfun(),
 * minqit(), minqir() in minqa/src/minqa.cpp). One MinqaCalfun lives on
 * the stack of each *_solve_c() call, so concurrent solves share nothing.
 *
 * Aborts (objfun error, non-finite x) unwind the Fortran-style goto code
 * via the MinqaAbort exception, caught in the extern "C" entry point. */

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <string>
#include "minqa_types.h"

struct MinqaAbort {
  int code;  // MINQA_ERR_OBJFUN or MINQA_ERR_NONFINITE_X
};

struct MinqaCalfun {
  minqa_c_objfun_t fn;
  void *userdata;
  int n;
  int iprint;
  minqa_print_t print;
  void *print_userdata;
  int feval = 0;

  MinqaCalfun(minqa_c_objfun_t fn_, void *userdata_, int n_,
              const minqa_options_t *opts)
      : fn(fn_), userdata(userdata_), n(n_), iprint(opts->iprint),
        print(opts->print), print_userdata(opts->print_userdata) {}

  void emit(const std::string &line) const {
    if (print != nullptr) print(line.c_str(), print_userdata);
  }

  static std::string fmt(const char *format, double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), format, v);
    return buf;
  }
  static std::string fmti(const char *format, int v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), format, v);
    return buf;
  }

  // F77_NAME(calfun): counts, rejects non-finite x, clamps non-finite f
  // to DBL_MAX, and prints per iprint exactly as minqa does.
  double operator()(const double *x) {
    feval++;
    for (int i = 0; i < n; i++)
      if (!std::isfinite(x[i])) throw MinqaAbort{MINQA_ERR_NONFINITE_X};
    double f = 0.0;
    if (fn(n, x, &f, userdata) < 0) throw MinqaAbort{MINQA_ERR_OBJFUN};
    if (!std::isfinite(f)) f = DBL_MAX;
    if (print != nullptr &&
        (iprint == 3 || (iprint > 3 && feval % iprint == 0))) {
      std::string line = fmti("%3d:", feval) + fmt("%#14.8g:", f);
      for (int i = 0; i < n; i++) line += fmt(" %#8g", x[i]);
      emit(line + "\n");
    }
    return f;
  }

  // F77_NAME(minqit): called when rho changes and iprint >= 2.
  // xbase and xopt are 0-based pointers to length-n arrays.
  void minqit(double rho, int nf, double fopt, const double *xbase,
              const double *xopt) const {
    if (print == nullptr || iprint < 2) return;
    std::string line = fmt("rho: %#8.2g", rho) + fmti(" eval: %3d", nf) +
                       fmt(" fn: %#12g par:", fopt);
    for (int i = 0; i < n; i++) line += fmt("%#8g ", xbase[i] + xopt[i]);
    emit(line + "\n");
  }

  // F77_NAME(minqir): output at return when iprint > 0.
  void minqir(double f, int nf, const double *x) const {
    if (print == nullptr || iprint <= 0) return;
    emit("At return\n");
    std::string line = fmti("eval: %3d", nf) + fmt(" fn: %#14.8g par:", f);
    for (int i = 0; i < n; i++) line += fmt(" %#8g", x[i]);
    emit(line + "\n");
  }
};

/* Maps the raw Fortran IERR to minqa.R's remapped code. */
inline int minqa_map_ierr(int raw) {
  switch (raw) {
    case 0: return MINQA_OK;
    case 390: return MINQA_ERR_MAXFUN;
    case 10: return MINQA_ERR_NPT;
    case 430: case 3701: return MINQA_ERR_NO_REDUCTION;
    case 20: return MINQA_ERR_RANGE;
    case 320: return MINQA_ERR_CANCELLATION;
    default: return raw;
  }
}

/* Shared argument checks for every *_solve_c() entry point. Returns
 * MINQA_OK, or zeroes *result, stores the code and returns it. */
inline int minqa_check_args(int n, const double *par, minqa_c_objfun_t objfun,
                            const minqa_options_t *opts,
                            minqa_result_t *result) {
  if (result == nullptr) return MINQA_ERR_INVALID;
  if (n < 0 || (n > 0 && par == nullptr) || objfun == nullptr ||
      !(opts->rhobeg > 0.0) || !(opts->rhoend > 0.0)) {
    minqa_result_zero(result);
    result->n = n;
    result->ierr = MINQA_ERR_INVALID;
    return MINQA_ERR_INVALID;
  }
  return MINQA_OK;
}

/* Common tail of every *_solve_c(): copy par, re-evaluate fval at par
 * (uncounted, unclamped, like minqa's rval()), map ierr. Returns
 * result->ierr. Defined in minqa_common.cpp. */
int minqa_finish(int n, const double *x, int raw_ierr, MinqaCalfun &calfun,
                 minqa_result_t *result);

/* Wraps a solver body with the shared error translation. body() must
 * run the port and return the raw Fortran IERR, leaving the final X in
 * x. Defined as a template so each solver's core stays in its own file. */
template <class Body>
int minqa_run(int n, double *x, MinqaCalfun &calfun, minqa_result_t *result,
              Body body) {
  minqa_result_zero(result);
  result->n = n;
  try {
    int raw = body();
    return minqa_finish(n, x, raw, calfun, result);
  } catch (const MinqaAbort &a) {
    result->ierr = a.code;
    result->feval = calfun.feval;
    return a.code;
  } catch (const std::bad_alloc &) {
    minqa_result_free(result);
    result->n = n;
    result->ierr = MINQA_ERR_NOMEM;
    return MINQA_ERR_NOMEM;
  } catch (...) {
    result->ierr = MINQA_ERR_INTERNAL;
    result->feval = calfun.feval;
    return MINQA_ERR_INTERNAL;
  }
}

#endif /* MINQA_COMMON_H */
