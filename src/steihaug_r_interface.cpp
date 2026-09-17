// R-objfun entry point for the Steihaug trust-region port, called by the
// shim in R/steihaug.R. Drives the same thread-safe steihaug_solve_c() a
// C/C++ consumer calls; only the callbacks (which evaluate R closures, so
// are main-thread only) differ.
#include <Rcpp.h>
#include <exception>
#include <string>
#include "steihaug_types.h"
using namespace Rcpp;

namespace {

struct SteihaugRObjfun {
  SEXP fn, gr, hess, hessvec;  // hess/hessvec may be R_NilValue
  std::exception_ptr error;
};

void steihaug_r_copy(SEXP v, int len, double *out, const char *what) {
  if (!Rf_isNumeric(v) || Rf_length(v) != len)
    Rcpp::stop("'%s' must return a numeric of length %d", what, len);
  NumericVector nv(v);
  std::copy(nv.begin(), nv.end(), out);
}

int steihaug_r_objfun(int n, const double *x, double *value, double *gradient,
                      double *hessian, void *userdata) {
  SteihaugRObjfun *d = static_cast<SteihaugRObjfun *>(userdata);
  try {
    NumericVector xx(x, x + n);
    if (value != nullptr) {
      Function fn(d->fn);
      SEXP v = fn(xx);
      if (!Rf_isNumeric(v) || Rf_length(v) != 1)
        Rcpp::stop("'fn' must return a single numeric value");
      *value = Rf_asReal(v);
    }
    if (gradient != nullptr) {
      Function gr(d->gr);
      steihaug_r_copy(gr(xx), n, gradient, "gr");
    }
    if (hessian != nullptr) {
      Function hs(d->hess);
      SEXP h = hs(xx);
      if (!Rf_isMatrix(h) || Rf_nrows(h) != n || Rf_ncols(h) != n)
        Rcpp::stop("'hess' must return a %d by %d numeric matrix", n, n);
      NumericMatrix hm(h);
      // The C contract is row-major.
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) hessian[i * n + j] = hm(i, j);
    }
    return 0;
  } catch (...) {
    d->error = std::current_exception();
    return -1;
  }
}

int steihaug_r_hessvec(int n, const double *x, const double *v, double *hv,
                       void *userdata) {
  SteihaugRObjfun *d = static_cast<SteihaugRObjfun *>(userdata);
  try {
    NumericVector xx(x, x + n), vv(v, v + n);
    Function hvf(d->hessvec);
    steihaug_r_copy(hvf(xx, vv), n, hv, "hessvec");
    return 0;
  } catch (...) {
    d->error = std::current_exception();
    return -1;
  }
}

const char *steihaug_status_name(int status) {
  switch (status) {
    case STEIHAUG_MAX_ITER: return "MaxIter";
    case STEIHAUG_MAX_COST_EVALS: return "MaxCostEvals";
    case STEIHAUG_MAX_GRADIENT_EVALS: return "MaxGradientEvals";
    case STEIHAUG_TARGET_COST: return "TargetCost";
    case STEIHAUG_GRADIENT_TOLERANCE: return "GradientTolerance";
    case STEIHAUG_RELATIVE_GRADIENT_TOLERANCE: return "RelativeGradientTolerance";
    case STEIHAUG_PARAM_TOLERANCE: return "ParamTolerance";
    case STEIHAUG_RELATIVE_PARAM_TOLERANCE: return "RelativeParamTolerance";
    case STEIHAUG_COST_TOLERANCE: return "CostTolerance";
    case STEIHAUG_RELATIVE_COST_TOLERANCE: return "RelativeCostTolerance";
    case STEIHAUG_SOLVER_CONVERGED: return "SolverConverged";
    case STEIHAUG_SOLVER_FAILED: return "SolverFailed";
    default: return "Error";
  }
}

double ctrl_dbl(List ctrl, const char *nm) { return as<double>(ctrl[nm]); }

// NULL in the R control list disables a tolerance / budget.
double ctrl_tol(List ctrl, const char *nm) {
  SEXP v = ctrl[nm];
  return Rf_isNull(v) ? STEIHAUG_DISABLED : as<double>(v);
}

}  // namespace

// [[Rcpp::export]]
List steihaug_solve_r(NumericVector par, Function fn, Function gr,
                      RObject hess, RObject hessvec, List control) {
  int n = par.size();
  SteihaugRObjfun d{fn, gr, hess, hessvec, nullptr};

  steihaug_options_t o = steihaug_options_default();
  o.initial_radius = ctrl_dbl(control, "rinit");
  o.max_radius = ctrl_dbl(control, "rmax");
  o.eta = ctrl_dbl(control, "eta");
  o.max_inner = as<int>(control["maxInner"]);
  o.cg_max_iter = as<int>(control["cgMaxit"]);
  o.kappa = ctrl_dbl(control, "kappa");
  o.theta = ctrl_dbl(control, "theta");
  o.max_iter = as<int>(control["maxit"]);
  SEXP mce = control["maxCostEvals"], mge = control["maxGradEvals"],
       tc = control["targetCost"];
  o.max_cost_evals = Rf_isNull(mce) ? -1 : static_cast<long long>(as<double>(mce));
  o.max_gradient_evals =
      Rf_isNull(mge) ? -1 : static_cast<long long>(as<double>(mge));
  o.has_target_cost = Rf_isNull(tc) ? 0 : 1;
  o.target_cost = Rf_isNull(tc) ? 0.0 : as<double>(tc);
  o.abs_gradient_tol = ctrl_tol(control, "gradTol");
  o.rel_gradient_tol = ctrl_tol(control, "relGradTol");
  o.abs_step_tol = ctrl_tol(control, "stepTol");
  o.rel_step_tol = ctrl_tol(control, "relStepTol");
  o.abs_cost_change_tol = ctrl_tol(control, "costTol");
  o.rel_cost_change_tol = ctrl_tol(control, "relCostTol");
  o.trace = as<bool>(control["trace"]) ? 1 : 0;

  steihaug_result_t res;
  steihaug_solve_c(n, par.begin(), steihaug_r_objfun,
                   hessvec.isNULL() ? nullptr : steihaug_r_hessvec, &d, &o,
                   &res);

  if (d.error) {
    steihaug_result_free(&res);
    std::rethrow_exception(d.error);
  }
  if (res.status == STEIHAUG_ERR_INVALID_OPTIONS) {
    steihaug_result_free(&res);
    Rcpp::stop("invalid steihaug() control settings");
  }
  if (res.status < 0 || res.argument == nullptr) {
    int code = res.status;
    steihaug_result_free(&res);
    Rcpp::stop("steihaug solver failed with internal error code %d", code);
  }

  NumericVector counts = NumericVector::create(
      _["fn"] = static_cast<double>(res.cost_evals),
      _["gr"] = static_cast<double>(res.gradient_evals),
      _["hess"] = static_cast<double>(res.hessian_evals),
      _["hessvec"] = static_cast<double>(res.hessvec_evals));
  List out = List::create(
      _["par"] = NumericVector(res.argument, res.argument + n),
      _["value"] = res.value,
      _["gradient"] = NumericVector(res.gradient, res.gradient + n),
      _["iterations"] = res.iterations, _["converged"] = res.converged != 0,
      _["status"] = res.status, _["message"] = steihaug_status_name(res.status),
      _["radius"] = res.radius, _["counts"] = counts);
  if (o.trace) {
    int m = res.trace_len;
    DataFrame tr = DataFrame::create(
        _["iter"] = IntegerVector(res.trace_iter, res.trace_iter + m),
        _["radius"] = NumericVector(res.trace_radius, res.trace_radius + m),
        _["value"] = NumericVector(res.trace_cost, res.trace_cost + m),
        _["gradnorm"] = NumericVector(res.trace_gradnorm, res.trace_gradnorm + m),
        _["trialValue"] =
            NumericVector(res.trace_trial_cost, res.trace_trial_cost + m),
        _["preddiff"] = NumericVector(res.trace_pred, res.trace_pred + m),
        _["rho"] = NumericVector(res.trace_rho, res.trace_rho + m),
        _["stepnorm"] = NumericVector(res.trace_step_norm, res.trace_step_norm + m),
        _["accept"] = LogicalVector(res.trace_accept, res.trace_accept + m),
        _["hitBoundary"] =
            LogicalVector(res.trace_hit_boundary, res.trace_hit_boundary + m),
        _["cgIter"] = IntegerVector(res.trace_cg_iter, res.trace_cg_iter + m));
    out["trace"] = tr;
  }
  steihaug_result_free(&res);
  return out;
}
