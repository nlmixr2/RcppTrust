#include "trust_types.h"
#include <cstdlib>
#include <cstring>

void trust_result_zero(trust_result_t *res) { std::memset(res, 0, sizeof(*res)); }

void trust_result_free(trust_result_t *res) {
  free(res->gradient);
  free(res->hessian);
  free(res->argument);
  free(res->argpath);
  free(res->argtry);
  free(res->steptype);
  free(res->accept);
  free(res->r);
  free(res->rho);
  free(res->valpath);
  free(res->valtry);
  free(res->preddiff);
  free(res->stepnorm);
  trust_result_zero(res);
}
