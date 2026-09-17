// R-objfun entry points for the BOBYQA/NEWUOA ports, called by the
// shims in R/minqa.R. These replace minqa's bobyqa_cpp/newuoa_cpp +
// calfun()/rval() glue, but drive the same thread-safe
// *_solve_c() cores a C/C++ consumer calls -- only the objective callback
// (which evaluates an R closure, so is main-thread only) differs.
#include <Rcpp.h>
#include <exception>
#include "minqa_types.h"
using namespace Rcpp;

namespace {

struct MinqaRObjfun {
  SEXP fn;
  std::exception_ptr error;
};

// Mirrors minqa's calfun(): Rf_asReal(fn(x)), evaluated in fn's
// environment. Any R error (or interrupt) is captured and re-thrown by the
// caller once the solver has unwound, so no longjmp crosses the core.
int minqa_r_objfun(int n, const double *x, double *f, void *userdata) {
  MinqaRObjfun *d = static_cast<MinqaRObjfun *>(userdata);
  try {
    NumericVector xx(x, x + n);
    Function fn(d->fn);
    SEXP v = fn(xx);
    // Coerce odd returns through R (unwind-protected by Rcpp) so a coercion
    // warning promoted to an error can't longjmp across the solver frames;
    // the value is the same NA-or-number Rf_asReal() gives minqa.
    if (!Rf_isNumeric(v) && !Rf_isLogical(v)) v = Function("as.double")(v);
    *f = Rf_asReal(v);
    return 0;
  } catch (...) {
    d->error = std::current_exception();
    return -1;
  }
}

void minqa_r_print(const char *line, void *) { Rprintf("%s", line); }

minqa_options_t minqa_r_options(Environment ctrl) {
  minqa_options_t opts;
  opts.npt = as<int>(ctrl.get("npt"));
  opts.rhobeg = as<double>(ctrl.get("rhobeg"));
  opts.rhoend = as<double>(ctrl.get("rhoend"));
  opts.iprint = as<int>(ctrl.get("iprint"));
  opts.maxfun = as<int>(ctrl.get("maxfun"));
  opts.print = minqa_r_print;
  opts.print_userdata = nullptr;
  return opts;
}

// Shape of minqa's rval() list (before the R shim adds `msg`).
List minqa_r_result(int n, minqa_result_t &res, MinqaRObjfun &d,
                    const char *cls) {
  if (d.error) {
    minqa_result_free(&res);
    std::rethrow_exception(d.error);
  }
  if (res.ierr == MINQA_ERR_NONFINITE_X) {
    minqa_result_free(&res);
    Rcpp::stop("non-finite x values not allowed in calfun");
  }
  if (res.ierr < 0 || res.par == nullptr) {
    int code = res.ierr;
    minqa_result_free(&res);
    Rcpp::stop("minqa solver failed with internal error code %d", code);
  }
  NumericVector par(res.par, res.par + n);
  double fval = res.fval;
  int raw = res.raw_ierr;
  IntegerVector feval = IntegerVector::create(res.feval);
  minqa_result_free(&res);
  List rr = List::create(_["par"] = par, _["fval"] = fval,
                         _["feval"] = feval, _["ierr"] = raw);
  rr.attr("class") = CharacterVector::create(cls, "minqa");
  return rr;
}

}  // namespace

// [[Rcpp::export]]
List minqa_bobyqa_r(NumericVector par, NumericVector lower, NumericVector upper,
                    Environment ctrl, Function fn) {
  int n = par.size();
  MinqaRObjfun d{fn, nullptr};
  minqa_options_t opts = minqa_r_options(ctrl);
  minqa_result_t res;
  bobyqa_solve_c(n, par.begin(), lower.begin(), upper.begin(), minqa_r_objfun,
                 &d, &opts, &res);
  return minqa_r_result(n, res, d, "bobyqa");
}

// [[Rcpp::export]]
List minqa_newuoa_r(NumericVector par, Environment ctrl, Function fn) {
  int n = par.size();
  MinqaRObjfun d{fn, nullptr};
  minqa_options_t opts = minqa_r_options(ctrl);
  minqa_result_t res;
  newuoa_solve_c(n, par.begin(), minqa_r_objfun, &d, &opts, &res);
  return minqa_r_result(n, res, d, "newuoa");
}
