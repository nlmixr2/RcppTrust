// Cross-package registration surface: builds the positionally-indexed
// list of function pointers consumed by inst/include/RcppTrust.h (the
// header-only pattern used throughout the nlmixr2 ecosystem, e.g.
// n1qn1c/inst/include/n1qn1c.h), and separately registers the same
// pointers via R_RegisterCCallable for completeness/introspection.
//
// Every solver AND its matching *_result_free must be exposed this way: a
// consumer never links against RcppTrust's shared library, so it cannot
// call either one directly by name -- only through the resolved pointer.
//
// Rcpp already auto-generates this package's R_init_RcppTrust() (in
// RcppExports.cpp, via useDynLib(.registration = TRUE)), so the
// R_RegisterCCallable() calls below are attached to that same
// initialization via `// [[Rcpp::init]]` rather than a hand-written
// src/init.c (which would collide with the generated one).
#include <Rcpp.h>
#include "minqa_types.h"
#include "steihaug_types.h"
#include "trust_types.h"
using namespace Rcpp;

// Slot order is part of the ABI consumed by inst/include/RcppTrust.h:
// only ever append, never reorder or remove.
// [[Rcpp::export]]
SEXP trust_ptr() {
  static const char *names[] = {"trust_solve_c",      "trust_result_free",
                                "bobyqa_solve_c",     "newuoa_solve_c",
                                "minqa_result_free",  "steihaug_solve_c",
                                "steihaug_result_free"};
  DL_FUNC fns[] = {(DL_FUNC)&trust_solve_c,     (DL_FUNC)&trust_result_free,
                   (DL_FUNC)&bobyqa_solve_c,    (DL_FUNC)&newuoa_solve_c,
                   (DL_FUNC)&minqa_result_free, (DL_FUNC)&steihaug_solve_c,
                   (DL_FUNC)&steihaug_result_free};
  const int nVec = sizeof(fns) / sizeof(fns[0]);
  SEXP ret = PROTECT(Rf_allocVector(VECSXP, nVec));
  SEXP retN = PROTECT(Rf_allocVector(STRSXP, nVec));
  for (int i = 0; i < nVec; i++) {
    SET_VECTOR_ELT(ret, i, R_MakeExternalPtrFn(fns[i], R_NilValue, R_NilValue));
    SET_STRING_ELT(retN, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(ret, R_NamesSymbol, retN);
  UNPROTECT(2);
  return ret;
}

// [[Rcpp::init]]
void trust_register_ccallable(DllInfo *dll) {
  R_RegisterCCallable("RcppTrust", "trust_solve_c", (DL_FUNC)&trust_solve_c);
  R_RegisterCCallable("RcppTrust", "trust_result_free", (DL_FUNC)&trust_result_free);
  R_RegisterCCallable("RcppTrust", "bobyqa_solve_c", (DL_FUNC)&bobyqa_solve_c);
  R_RegisterCCallable("RcppTrust", "newuoa_solve_c", (DL_FUNC)&newuoa_solve_c);
  R_RegisterCCallable("RcppTrust", "minqa_result_free", (DL_FUNC)&minqa_result_free);
  R_RegisterCCallable("RcppTrust", "steihaug_solve_c", (DL_FUNC)&steihaug_solve_c);
  R_RegisterCCallable("RcppTrust", "steihaug_result_free",
                      (DL_FUNC)&steihaug_result_free);
}
