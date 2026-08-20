#ifndef TRUST_ROBJFUN_H
#define TRUST_ROBJFUN_H

#include <RcppArmadillo.h>
#include "trust_core.h"

/* Evaluator for the R-objfun path: each call goes through the tiny R
 * helper `.trustTryEval <- function(f, theta) try(f(theta), silent =
 * TRUE)`, so that an error raised inside the user's R objfun produces a
 * genuine R "try-error" object (with its "condition" attribute) exactly
 * as upstream trust()'s own try(objfun(theta, ...)) does -- this is the
 * one piece of real control flow kept in R, needed only because R
 * condition objects can't be reconstructed after crossing into C++.
 *
 * Not thread-safe (touches the R API via Rcpp::Function on every call);
 * only ever instantiated on the main R thread. */
struct RObjfunEvaluator {
  Rcpp::Function tryEval;  // .trustTryEval
  Rcpp::Function objfun;   // already wrapped to be unary: function(theta)
  int dimen;
  bool minimize;

  Rcpp::RObject lastErrorObject;

  RObjfunEvaluator(Rcpp::Function tryEval_, Rcpp::Function objfun_, int dimen_,
                    bool minimize_)
      : tryEval(tryEval_), objfun(objfun_), dimen(dimen_), minimize(minimize_) {}

  TrustEvalStatus operator()(const arma::vec &theta, TrustEvalOutput &out);
};

#endif /* TRUST_ROBJFUN_H */
