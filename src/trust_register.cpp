// Cross-package registration surface: builds the positionally-indexed
// list of function pointers consumed by inst/include/RcppTrust.h (the
// header-only pattern used throughout the nlmixr2 ecosystem, e.g.
// n1qn1c/inst/include/n1qn1c.h), and separately registers the same
// pointers via R_RegisterCCallable for completeness/introspection.
//
// Both trust_solve_c AND trust_result_free must be exposed this way: a
// consumer never links against RcppTrust's shared library, so it cannot
// call either one directly by name -- only through the resolved pointer.
//
// Rcpp already auto-generates this package's R_init_RcppTrust() (in
// RcppExports.cpp, via useDynLib(.registration = TRUE)), so the
// R_RegisterCCallable() calls below are attached to that same
// initialization via `// [[Rcpp::init]]` rather than a hand-written
// src/init.c (which would collide with the generated one).
#include <Rcpp.h>
#include "trust_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP trust_ptr() {
#define nVec 2
  SEXP p0 = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&trust_solve_c, R_NilValue, R_NilValue));
  SEXP p1 =
      PROTECT(R_MakeExternalPtrFn((DL_FUNC)&trust_result_free, R_NilValue, R_NilValue));
  SEXP ret = PROTECT(Rf_allocVector(VECSXP, nVec));
  SEXP retN = PROTECT(Rf_allocVector(STRSXP, nVec));
  SET_VECTOR_ELT(ret, 0, p0);
  SET_STRING_ELT(retN, 0, Rf_mkChar("trust_solve_c"));
  SET_VECTOR_ELT(ret, 1, p1);
  SET_STRING_ELT(retN, 1, Rf_mkChar("trust_result_free"));
  Rf_setAttrib(ret, R_NamesSymbol, retN);
  UNPROTECT(nVec + 2);
#undef nVec
  return ret;
}

// [[Rcpp::init]]
void trust_register_ccallable(DllInfo *dll) {
  R_RegisterCCallable("RcppTrust", "trust_solve_c", (DL_FUNC)&trust_solve_c);
  R_RegisterCCallable("RcppTrust", "trust_result_free", (DL_FUNC)&trust_result_free);
}
