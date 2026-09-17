#include "minqa_common.h"
#include <cstdlib>
#include <cstring>
#include <new>

void minqa_result_zero(minqa_result_t *res) { std::memset(res, 0, sizeof(*res)); }

void minqa_result_free(minqa_result_t *res) {
  free(res->par);
  minqa_result_zero(res);
}

int minqa_finish(int n, const double *x, int raw_ierr, MinqaCalfun &calfun,
                 minqa_result_t *result) {
  result->raw_ierr = raw_ierr;
  result->ierr = minqa_map_ierr(raw_ierr);
  result->feval = calfun.feval;
  size_t bytes = n > 0 ? static_cast<size_t>(n) * sizeof(double) : sizeof(double);
  result->par = static_cast<double *>(malloc(bytes));
  if (result->par == nullptr) throw std::bad_alloc();
  if (n > 0) std::memcpy(result->par, x, static_cast<size_t>(n) * sizeof(double));
  double f = 0.0;
  if (calfun.fn(n, x, &f, calfun.userdata) < 0) {
    result->ierr = MINQA_ERR_OBJFUN;
    return result->ierr;
  }
  result->fval = f;
  return result->ierr;
}
