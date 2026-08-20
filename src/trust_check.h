#ifndef TRUST_CHECK_H
#define TRUST_CHECK_H

#include <string>

/* Numeric-only port of upstream trust's check.objfun.output(), restricted
 * to the checks that make sense once value/gradient/hessian are already
 * plain C arrays (i.e. everything except the R-type/shape checks, which
 * only apply to the R-objfun path and are done directly on the R list in
 * trust_robjfun.cpp). Used by both the R and C objfun paths so the
 * "is this evaluation acceptable" logic is identical either way.
 *
 * gradient/hessian may be NULL only when value is not finite (mirrors
 * check.objfun.output's `if (is.finite(foo)) { check gradient/hessian }`).
 *
 * Returns true if ok; on failure returns false and sets errmsg to the
 * exact upstream stop() message text. */
bool trust_check_numeric_output(double value, const double *gradient,
                                 const double *hessian, int dimen,
                                 bool minimize, std::string &errmsg);

#endif /* TRUST_CHECK_H */
