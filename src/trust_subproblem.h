#ifndef TRUST_SUBPROBLEM_H
#define TRUST_SUBPROBLEM_H

#include <RcppArmadillo.h>
#include "trust_types.h"

/* Direct port of the trust-region subproblem solve in upstream trust()'s
 * main loop (the block from "solve trust region subproblem" through the
 * hard-hard case), operating on an already-computed symmetric
 * eigendecomposition of B (eigenvalues ascending or descending, order
 * does not matter to this algorithm) and the gradient rotated into the
 * eigenbasis, gq = t(eigenvectors) %*% g.
 *
 * eigval, eigvec: eigendecomposition of the (possibly rescaled/negated)
 *   Hessian B, i.e. B == eigvec * diag(eigval) * eigvec.t().
 * gq: eigenvectors.t() * g, length d.
 * r: current trust region radius.
 *
 * Returns ptry (the trial step, length d) and, via the out-params,
 * which case was taken -- this mirrors is.newton/is.hard/is.easy and the
 * derived steptype exactly. */
struct TrustSubproblemResult {
  arma::vec ptry;
  bool is_newton;
  bool is_hard;
  bool is_easy;
  trust_steptype_t steptype;
};

TrustSubproblemResult trust_solve_subproblem(const arma::vec &eigval,
                                              const arma::mat &eigvec,
                                              const arma::vec &gq, double r);

#endif /* TRUST_SUBPROBLEM_H */
