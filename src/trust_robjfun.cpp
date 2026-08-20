#include "trust_robjfun.h"
#include <cmath>
#include <string>
using namespace Rcpp;

namespace {

// Direct port of upstream check.objfun.output(), operating on the R list
// objfun returned (rather than on plain doubles like
// trust_check_numeric_output(), which this intentionally duplicates for
// the R-shaped checks -- the R path needs the exact type/shape checks
// [is.list/is.numeric/is.matrix/length/dim] that only make sense on R
// objects, in upstream's exact order, producing upstream's exact
// messages). Throws TrustValidationError on any failure, matching
// check.objfun.output()'s uncaught stop() semantics.
void check_objfun_list(SEXP obj, int dimen, bool minimize, double &value,
                        arma::vec &gradient, arma::mat &hessian,
                        bool &value_finite) {
  if (!Rf_isNewList(obj))
    throw TrustValidationError("objfun returned object that is not a list");
  List lst(obj);

  if (!lst.containsElementNamed("value"))
    throw TrustValidationError(
        "objfun returned list that does not have a component 'value'");
  SEXP vsxp = lst["value"];
  if (!Rf_isNumeric(vsxp))
    throw TrustValidationError("objfun returned value that is not numeric");
  if (Rf_length(vsxp) != 1)
    throw TrustValidationError("objfun returned value that is not scalar");
  double val = as<double>(vsxp);
  if (ISNAN(val))  // R's is.na()||is.nan() for a double: both covered by ISNAN
    throw TrustValidationError("objfun returned value that is NA or NaN");
  if (minimize && val == R_NegInf)
    throw TrustValidationError("objfun returned -Inf value in minimization");
  if ((!minimize) && val == R_PosInf)
    throw TrustValidationError("objfun returned +Inf value in maximization");

  value = val;
  value_finite = R_FINITE(val) != 0;

  if (value_finite) {
    if (!lst.containsElementNamed("gradient"))
      throw TrustValidationError(
          "objfun returned list without component 'gradient' when value is "
          "finite");
    SEXP gsxp = lst["gradient"];
    if (!Rf_isNumeric(gsxp))
      throw TrustValidationError("objfun returned gradient that is not numeric");
    if (Rf_length(gsxp) != dimen)
      throw TrustValidationError(
          "objfun returned gradient that is not vector of length " +
          std::to_string(dimen));
    NumericVector gvec(gsxp);
    for (double gi : gvec)
      if (!R_FINITE(gi))
        throw TrustValidationError(
            "objfun returned gradient not having all elements finite");
    gradient = arma::vec(gvec.begin(), gvec.size());

    if (!lst.containsElementNamed("hessian"))
      throw TrustValidationError(
          "objfun returned list without component 'hessian' when value is "
          "finite");
    SEXP hsxp = lst["hessian"];
    if (!Rf_isNumeric(hsxp))
      throw TrustValidationError("objfun returned hessian that is not numeric");
    if (!Rf_isMatrix(hsxp))
      throw TrustValidationError("objfun returned hessian that is not matrix");
    NumericMatrix hmat(hsxp);
    if (hmat.nrow() != dimen || hmat.ncol() != dimen)
      throw TrustValidationError(
          "objfun returned hessian that is not " + std::to_string(dimen) +
          " by " + std::to_string(dimen) + " matrix");
    for (double hi : hmat)
      if (!R_FINITE(hi))
        throw TrustValidationError(
            "objfun returned hessian not having all elements finite");
    hessian = arma::mat(hmat.begin(), dimen, dimen);
  }
}

}  // namespace

TrustEvalStatus RObjfunEvaluator::operator()(const arma::vec &theta,
                                              TrustEvalOutput &out) {
  NumericVector th(theta.begin(), theta.end());
  RObject result = tryEval(objfun, th);
  if (Rf_inherits(result, "try-error")) {
    lastErrorObject = result;
    return TrustEvalStatus::ERROR;
  }
  double value;
  bool value_finite;
  check_objfun_list(result, dimen, minimize, value, out.gradient, out.hessian,
                     value_finite);
  out.value = value;
  return value_finite ? TrustEvalStatus::OK : TrustEvalStatus::INFEASIBLE;
}
