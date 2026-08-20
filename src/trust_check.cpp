#include "trust_check.h"
#include <cmath>
#include <sstream>

bool trust_check_numeric_output(double value, const double *gradient,
                                 const double *hessian, int dimen,
                                 bool minimize, std::string &errmsg) {
  if (std::isnan(value)) {
    // R's is.na() is true for NaN too, and is.na() is checked first.
    errmsg = "objfun returned value that is NA or NaN";
    return false;
  }
  if (minimize && value == -HUGE_VAL) {
    errmsg = "objfun returned -Inf value in minimization";
    return false;
  }
  if ((!minimize) && value == HUGE_VAL) {
    errmsg = "objfun returned +Inf value in maximization";
    return false;
  }
  if (std::isfinite(value)) {
    if (gradient == nullptr) {
      errmsg =
          "objfun returned list without component 'gradient' when value is "
          "finite";
      return false;
    }
    for (int i = 0; i < dimen; i++) {
      if (!std::isfinite(gradient[i])) {
        errmsg = "objfun returned gradient not having all elements finite";
        return false;
      }
    }
    if (hessian == nullptr) {
      errmsg =
          "objfun returned list without component 'hessian' when value is "
          "finite";
      return false;
    }
    for (int i = 0; i < dimen * dimen; i++) {
      if (!std::isfinite(hessian[i])) {
        errmsg = "objfun returned hessian not having all elements finite";
        return false;
      }
    }
  }
  return true;
}
