#include "trust_core.h"
#include "trust_types.h"
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

// Thread-safe Evaluator wrapping the C objfun contract: no R API calls,
// safe to instantiate/run concurrently from independent threads as long
// as each thread uses its own CObjfunEvaluator/userdata.
struct CObjfunEvaluator {
  trust_c_objfun_t fn;
  void *userdata;
  int n;

  TrustEvalStatus operator()(const arma::vec &theta, TrustEvalOutput &out) {
    std::vector<double> grad(n), hess(static_cast<size_t>(n) * n);
    double value = 0.0;
    int rc = fn(n, theta.memptr(), &value, grad.data(), hess.data(), userdata);
    if (rc < 0) return TrustEvalStatus::ERROR;
    if (rc == 1) {
      out.value = value;
      out.gradient.reset();
      out.hessian.reset();
      return TrustEvalStatus::INFEASIBLE;
    }
    out.value = value;
    out.gradient = arma::vec(grad.data(), n);
    // Symmetric by construction (Hessian), so row-major vs column-major
    // storage of the caller's buffer is immaterial here.
    out.hessian = arma::mat(hess.data(), n, n);
    return TrustEvalStatus::OK;
  }
};

double *dup_doubles(const double *src, size_t count) {
  size_t bytes = count > 0 ? count * sizeof(double) : sizeof(double);
  double *p = static_cast<double *>(malloc(bytes));
  if (p == nullptr) throw std::bad_alloc();
  if (src != nullptr && count > 0) memcpy(p, src, count * sizeof(double));
  return p;
}

double *dup_vec(const arma::vec &v) { return dup_doubles(v.memptr(), v.n_elem); }

// Row-major flatten of a list of arma::vec (one row per entry) into a
// freshly malloc'd buffer, matching trust_result_t's documented layout.
double *dup_rows(const std::vector<arma::vec> &rows, int n) {
  size_t count = rows.size() * static_cast<size_t>(n);
  size_t bytes = count > 0 ? count * sizeof(double) : sizeof(double);
  double *p = static_cast<double *>(malloc(bytes));
  if (p == nullptr) throw std::bad_alloc();
  for (size_t i = 0; i < rows.size(); i++)
    for (int j = 0; j < n; j++) p[i * n + j] = rows[i](j);
  return p;
}

int *dup_ints(const std::vector<int> &src) {
  size_t count = src.size();
  size_t bytes = count > 0 ? count * sizeof(int) : sizeof(int);
  int *p = static_cast<int *>(malloc(bytes));
  if (p == nullptr) throw std::bad_alloc();
  if (count > 0) memcpy(p, src.data(), count * sizeof(int));
  return p;
}

int *dup_steptypes(const std::vector<trust_steptype_t> &src) {
  size_t count = src.size();
  size_t bytes = count > 0 ? count * sizeof(int) : sizeof(int);
  int *p = static_cast<int *>(malloc(bytes));
  if (p == nullptr) throw std::bad_alloc();
  for (size_t i = 0; i < count; i++) p[i] = static_cast<int>(src[i]);
  return p;
}


}  // namespace

extern "C" int trust_solve_c(int n, const double *parinit, trust_c_objfun_t objfun,
                              void *userdata, const trust_options_t *opts,
                              trust_result_t *result) {
  trust_result_zero(result);
  result->n = n;

  try {
    arma::vec p0(parinit, n);
    CObjfunEvaluator eval{objfun, userdata, n};

    arma::vec parscale_vec;
    const arma::vec *parscale_ptr = nullptr;
    if (opts->has_parscale) {
      parscale_vec = arma::vec(opts->parscale, n);
      parscale_ptr = &parscale_vec;
    }

    TrustRunResult run =
        trust_core_run(eval, p0, opts->rinit, opts->rmax, parscale_ptr,
                        opts->iterlim, opts->fterm, opts->mterm,
                        opts->minimize != 0, opts->blather != 0);

    result->iterations = run.iterations;

    if (run.error_occurred) {
      result->error = run.error_source;  // 1 = first call, 2 = mid-loop, 3 = final call
      result->converged = 0;
      if (run.error_source == 2) {
        result->argument = dup_vec(run.argument_at_error);
      } else if (run.error_source == 3) {
        result->argument = dup_vec(run.argument);
        result->converged = run.converged ? 1 : 0;
      }
      // error_source == 1: no argument recorded beyond parinit (caller
      // already has it), matching upstream's early-return shape.
      return result->error;
    }

    result->converged = run.converged ? 1 : 0;
    result->argument = dup_vec(run.argument);
    if (run.have_final_eval) {
      result->value = run.value;
      result->gradient = dup_vec(run.gradient);
      result->hessian = dup_doubles(run.hessian.memptr(), run.hessian.n_elem);
    }

    if (opts->blather) {
      result->path_len = static_cast<int>(run.argpath.size());
      result->try_len = static_cast<int>(run.argtry.size());
      result->argpath = dup_rows(run.argpath, n);
      result->argtry = dup_rows(run.argtry, n);
      result->r = dup_doubles(run.r.data(), run.r.size());
      result->valpath = dup_doubles(run.valpath.data(), run.valpath.size());
      result->rho = dup_doubles(run.rho.data(), run.rho.size());
      result->valtry = dup_doubles(run.valtry.data(), run.valtry.size());
      result->preddiff = dup_doubles(run.preddiff.data(), run.preddiff.size());
      result->stepnorm = dup_doubles(run.stepnorm.data(), run.stepnorm.size());
      result->accept = dup_ints(run.accept);
      result->steptype = dup_steptypes(run.steptype);
      if (!opts->minimize) {
        for (int i = 0; i < result->try_len; i++) result->preddiff[i] = -result->preddiff[i];
      }
    }

    return 0;
  } catch (const TrustValidationError &) {
    result->error = -2;
    return -2;
  } catch (const TrustInfeasibleStartError &) {
    result->error = -3;
    return -3;
  } catch (const std::bad_alloc &) {
    trust_result_free(result);
    result->error = -4;
    return -4;
  } catch (...) {
    result->error = -99;
    return -99;
  }
}
