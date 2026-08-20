// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "trust_check.h"
#include "trust_core.h"
#include "trust_robjfun.h"
#include "trust_rootfind.h"
#include "trust_subproblem.h"
using namespace Rcpp;

namespace {

List blather_to_list(const TrustRunResult &run, int n) {
  int path_len = static_cast<int>(run.argpath.size());
  int try_len = static_cast<int>(run.argtry.size());
  NumericMatrix argpath(path_len, n), argtry(try_len, n);
  for (int i = 0; i < path_len; i++)
    for (int j = 0; j < n; j++) argpath(i, j) = run.argpath[i](j);
  for (int i = 0; i < try_len; i++)
    for (int j = 0; j < n; j++) argtry(i, j) = run.argtry[i](j);

  static const char *names[] = {"Newton", "easy-easy", "hard-easy", "hard-hard"};
  CharacterVector steptype(try_len);
  LogicalVector accept(try_len);
  for (int i = 0; i < try_len; i++) {
    steptype[i] = names[static_cast<int>(run.steptype[i])];
    accept[i] = run.accept[i] != 0;
  }

  return List::create(
      _["argpath"] = argpath, _["argtry"] = argtry, _["steptype"] = steptype,
      _["accept"] = accept, _["r"] = wrap(run.r), _["rho"] = wrap(run.rho),
      _["valpath"] = wrap(run.valpath), _["valtry"] = wrap(run.valtry),
      _["preddiff"] = wrap(run.preddiff), _["stepnorm"] = wrap(run.stepnorm));
}

}  // namespace

// Entry point called by R/trust.R's trust() shim. objfun1 must already be
// a unary function of theta (the R shim wraps any `...` extra arguments
// into it), and tryEval must be the internal .trustTryEval helper.
// [[Rcpp::export]]
List trust_solve_r(Function objfun1, Function tryEval, NumericVector parinit,
                    double rinit, double rmax, Nullable<NumericVector> parscale,
                    int iterlim, double fterm, double mterm, bool minimize,
                    bool blather) {
  int n = parinit.size();
  arma::vec p0(parinit.begin(), n);

  arma::vec parscale_vec;
  const arma::vec *parscale_ptr = nullptr;
  if (parscale.isNotNull()) {
    NumericVector ps(parscale);
    parscale_vec = arma::vec(ps.begin(), ps.size());
    parscale_ptr = &parscale_vec;
  }

  RObjfunEvaluator eval(tryEval, objfun1, n, minimize);
  TrustRunResult run = trust_core_run(eval, p0, rinit, rmax, parscale_ptr, iterlim,
                                       fterm, mterm, minimize, blather);

  if (run.error_occurred) {
    if (run.error_source == 1) {
      Rcpp::warning("error in first call to objfun");
      return List::create(_["error"] = eval.lastErrorObject, _["argument"] = parinit,
                           _["converged"] = false, _["iterations"] = 0);
    }
    List out;
    out["error"] = eval.lastErrorObject;
    if (run.error_source == 2) {
      out["argument"] = NumericVector(run.argument_at_error.begin(),
                                       run.argument_at_error.end());
      out["converged"] = false;
    } else {  // error_source == 3: error only on the final post-loop re-evaluation
      out["argument"] = NumericVector(run.argument.begin(), run.argument.end());
      out["converged"] = run.converged;
    }
    out["iterations"] = run.iterations;
    if (blather) {
      List bl = blather_to_list(run, n);
      CharacterVector blNames = bl.names();
      for (int i = 0; i < bl.size(); i++) out[as<std::string>(blNames[i])] = bl[i];
    }
    if (run.error_source == 3) Rcpp::warning("error in last call to objfun");
    return out;
  }

  List out;
  out["value"] = run.value;
  out["gradient"] = NumericVector(run.gradient.begin(), run.gradient.end());
  out["hessian"] = wrap(run.hessian);
  out["argument"] = NumericVector(run.argument.begin(), run.argument.end());
  out["converged"] = run.converged;
  out["iterations"] = run.iterations;
  if (blather) {
    List bl = blather_to_list(run, n);
    CharacterVector blNames = bl.names();
    for (int i = 0; i < bl.size(); i++) out[as<std::string>(blNames[i])] = bl[i];
  }
  return out;
}

// Thin test shim for trust_check_numeric_output(), used only by
// tests/testthat/test-subproblem-unit.R while building the port
// incrementally; exercises the exact numeric validation logic shared by
// both the R and C objfun paths.
// [[Rcpp::export]]
List trust_check_numeric_output_test(double value, Nullable<NumericVector> gradient,
                                      Nullable<NumericVector> hessian, int dimen,
                                      bool minimize) {
  std::string msg;
  const double *g = gradient.isNull() ? nullptr : NumericVector(gradient).begin();
  const double *h = hessian.isNull() ? nullptr : NumericVector(hessian).begin();
  bool ok = trust_check_numeric_output(value, g, h, dimen, minimize, msg);
  return List::create(_["ok"] = ok, _["message"] = msg);
}

// Thin test shim for trust_solve_lagrange_multiplier(), used only while
// building the port incrementally.
// [[Rcpp::export]]
double trust_solve_lagrange_multiplier_test(NumericVector beta, NumericVector gq,
                                             double r, double C1, double C2,
                                             double C3) {
  arma::vec b(beta.begin(), beta.size(), false);
  arma::vec g(gq.begin(), gq.size(), false);
  return trust_solve_lagrange_multiplier(b, g, r, C1, C2, C3);
}

// Thin test shim for trust_solve_subproblem(), used only while building
// the port incrementally.
// [[Rcpp::export]]
List trust_solve_subproblem_test(NumericVector eigval, NumericMatrix eigvec,
                                  NumericVector gq, double r) {
  arma::vec ev(eigval.begin(), eigval.size(), false);
  arma::mat V(eigvec.begin(), eigvec.nrow(), eigvec.ncol(), false);
  arma::vec g(gq.begin(), gq.size(), false);
  TrustSubproblemResult res = trust_solve_subproblem(ev, V, g, r);
  static const char *names[] = {"Newton", "easy-easy", "hard-easy", "hard-hard"};
  return List::create(
      _["ptry"] = NumericVector(res.ptry.begin(), res.ptry.end()),
      _["is_newton"] = res.is_newton, _["is_hard"] = res.is_hard,
      _["is_easy"] = res.is_easy, _["steptype"] = names[(int)res.steptype]);
}
