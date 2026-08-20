#ifndef TRUST_CORE_H
#define TRUST_CORE_H

#include <RcppArmadillo.h>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "trust_check.h"
#include "trust_subproblem.h"

/* Shared iteration loop, ported line-for-line from upstream trust()'s
 * main loop, parameterized over an Evaluator so the exact same control
 * flow and numeric algorithm runs whether the objective function is a
 * thread-safe C callback or an R closure. Only the Evaluator differs
 * between the two paths (see trust_types.h's CObjfunEvaluator in
 * trust_core.cpp, and trust_robjfun.h's RObjfunEvaluator).
 *
 * Evaluator must be callable as:
 *   TrustEvalStatus operator()(const arma::vec &theta, TrustEvalOutput &out)
 * and must itself remember any error detail needed later (an R condition
 * object for the R path, or nothing beyond a status code for the C path)
 * -- trust_core_run() only tracks *that* an error occurred and where. */

enum class TrustEvalStatus { OK, INFEASIBLE, ERROR };

struct TrustEvalOutput {
  double value = 0.0;
  arma::vec gradient;  // valid only when status == OK
  arma::mat hessian;   // valid only when status == OK
};

// Thrown when objfun() succeeds but its *return value* is invalid (the
// check.objfun.output() checks) -- upstream lets this propagate
// uncaught out of trust(), unlike an error raised inside objfun() itself
// (which try() catches gracefully). Callers translate this to whatever
// "hard, uncaught failure" mechanism fits their surface (Rcpp::stop() for
// the R path; a negative error code for the thread-safe C path).
struct TrustValidationError : std::runtime_error {
  explicit TrustValidationError(const std::string &msg) : std::runtime_error(msg) {}
};

// Thrown when parinit itself is infeasible (a value must be finite there,
// even though later trial points are allowed to signal infeasibility via
// +-Inf). Also an uncaught, non-graceful failure upstream.
struct TrustInfeasibleStartError : std::runtime_error {
  TrustInfeasibleStartError() : std::runtime_error("parinit not feasible") {}
};

struct TrustRunResult {
  bool error_occurred = false;
  int error_source = 0;  // 1 = first call, 2 = mid-loop (theta.try), 3 = final call
  bool converged = false;
  int iterations = 0;

  // Present when !error_occurred, or when error_source == 3 (error only
  // on the final re-evaluation -- argument/converged are still valid).
  arma::vec argument;
  double value = 0.0;
  arma::vec gradient;
  arma::mat hessian;
  bool have_final_eval = false;

  // Present when error_occurred && error_source == 2: the trial point
  // being evaluated when the error occurred (R's `theta.try`).
  arma::vec argument_at_error;

  bool blather = false;
  std::vector<arma::vec> argpath, argtry;  // path_len, try_len respectively
  std::vector<trust_steptype_t> steptype;  // try_len
  std::vector<int> accept;                 // try_len (0/1)
  std::vector<double> r, valpath;          // path_len
  std::vector<double> rho, valtry, preddiff, stepnorm;  // try_len
};

template <class Evaluator>
TrustRunResult trust_core_run(Evaluator &eval, const arma::vec &parinit,
                               double rinit, double rmax,
                               const arma::vec *parscale, int iterlim,
                               double fterm, double mterm, bool minimize,
                               bool blather) {
  const int d = static_cast<int>(parinit.n_elem);
  TrustRunResult res;
  res.blather = blather;

  double r = rinit;
  arma::vec theta = parinit;

  TrustEvalOutput out;
  TrustEvalStatus status = eval(theta, out);
  if (status == TrustEvalStatus::ERROR) {
    res.error_occurred = true;
    res.error_source = 1;
    res.argument_at_error = theta;
    res.iterations = 0;
    return res;
  }
  {
    std::string msg;
    const double *g = (status == TrustEvalStatus::OK) ? out.gradient.memptr() : nullptr;
    const double *h = (status == TrustEvalStatus::OK) ? out.hessian.memptr() : nullptr;
    if (!trust_check_numeric_output(out.value, g, h, d, minimize, msg))
      throw TrustValidationError(msg);
  }
  if (!std::isfinite(out.value)) throw TrustInfeasibleStartError();

  bool accept = true;
  double out_value = out.value;   // value at the most recent evaluation
  double value_save = out.value;  // value at the last *accepted* point

  arma::vec B_g;   // g, possibly rescaled/negated (persists across rejects)
  arma::mat B_h;   // B, possibly rescaled/negated
  double f = 0.0;
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec gq;

  bool error_break = false;
  bool terminated = false;
  arma::vec theta_try;
  TrustEvalOutput out_try;

  int iterations_ran = 0;
  for (int iiter = 1; iiter <= iterlim; ++iiter) {
    iterations_ran = iiter;

    if (blather) {
      res.argpath.push_back(theta);
      res.r.push_back(r);
      res.valpath.push_back(accept ? out_value : value_save);
    }

    if (accept) {
      B_h = out.hessian;
      B_g = out.gradient;
      f = out.value;
      value_save = f;
      if (parscale != nullptr) {
        B_h = B_h / ((*parscale) * (*parscale).t());
        B_g = B_g / (*parscale);
      }
      if (!minimize) {
        B_h = -B_h;
        B_g = -B_g;
        f = -f;
      }
      arma::eig_sym(eigval, eigvec, B_h);
      gq = eigvec.t() * B_g;
    }

    TrustSubproblemResult sub = trust_solve_subproblem(eigval, eigvec, gq, r);
    const arma::vec &ptry = sub.ptry;

    double preddiff = arma::dot(ptry, B_g + (B_h * ptry) / 2.0);
    if (parscale != nullptr) {
      theta_try = theta + ptry / (*parscale);
    } else {
      theta_try = theta + ptry;
    }

    TrustEvalStatus status_try = eval(theta_try, out_try);
    if (status_try == TrustEvalStatus::ERROR) {
      error_break = true;
      break;
    }
    {
      std::string msg;
      const double *g =
          (status_try == TrustEvalStatus::OK) ? out_try.gradient.memptr() : nullptr;
      const double *h =
          (status_try == TrustEvalStatus::OK) ? out_try.hessian.memptr() : nullptr;
      if (!trust_check_numeric_output(out_try.value, g, h, d, minimize, msg))
        throw TrustValidationError(msg);
    }
    out = out_try;
    out_value = out.value;

    double ftry = out.value;
    if (!minimize) ftry = -ftry;
    double rho;
    bool is_terminate;
    if (ftry < HUGE_VAL) {
      is_terminate = std::fabs(ftry - f) < fterm || std::fabs(preddiff) < mterm;
      rho = (ftry - f) / preddiff;
    } else {
      is_terminate = false;
      rho = -HUGE_VAL;
    }

    if (is_terminate) {
      if (ftry < f) {
        accept = true;
        theta = theta_try;
      }
      // else: accept stays whatever it already was from the previous
      // iteration's decision (matches R: no reassignment in this branch).
    } else {
      if (rho < 0.25) {
        accept = false;
        r = r / 4.0;
      } else {
        accept = true;
        theta = theta_try;
        if (rho > 0.75 && !sub.is_newton) r = std::min(2.0 * r, rmax);
      }
    }

    if (blather) {
      res.argtry.push_back(theta_try);
      res.valtry.push_back(out.value);
      res.accept.push_back(accept ? 1 : 0);
      res.preddiff.push_back(preddiff);
      res.stepnorm.push_back(arma::norm(ptry, 2));
      res.steptype.push_back(sub.steptype);
      res.rho.push_back(rho);
    }

    terminated = is_terminate;
    if (is_terminate) break;
  }
  res.iterations = iterations_ran;

  if (error_break) {
    res.error_occurred = true;
    res.error_source = 2;
    res.argument_at_error = theta_try;
    res.converged = false;
  } else {
    TrustEvalOutput out_final;
    TrustEvalStatus status_final = eval(theta, out_final);
    if (status_final == TrustEvalStatus::ERROR) {
      res.error_occurred = true;
      res.error_source = 3;
    } else {
      std::string msg;
      const double *g =
          (status_final == TrustEvalStatus::OK) ? out_final.gradient.memptr() : nullptr;
      const double *h =
          (status_final == TrustEvalStatus::OK) ? out_final.hessian.memptr() : nullptr;
      if (!trust_check_numeric_output(out_final.value, g, h, d, minimize, msg))
        throw TrustValidationError(msg);
      res.have_final_eval = true;
      res.value = out_final.value;
      res.gradient = out_final.gradient;
      res.hessian = out_final.hessian;
    }
    res.argument = theta;
    res.converged = terminated;
  }

  return res;
}

#endif /* TRUST_CORE_H */
