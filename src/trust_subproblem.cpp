#include "trust_subproblem.h"
#include "trust_rootfind.h"
#include <cmath>

TrustSubproblemResult trust_solve_subproblem(const arma::vec &eigval,
                                              const arma::mat &eigvec,
                                              const arma::vec &gq, double r) {
  TrustSubproblemResult res;
  bool is_newton = false;
  arma::vec ptry;

  // ##### try for Newton #####
  if (arma::all(eigval > 0.0)) {
    arma::vec wtry = gq / eigval;  // elementwise
    ptry = -(eigvec * wtry);
    if (arma::norm(ptry, 2) <= r) is_newton = true;
  }

  bool is_hard = false, is_easy = false;
  trust_steptype_t steptype = TRUST_STEP_NEWTON;

  // ##### non-Newton #####
  if (!is_newton) {
    double lambda_min = eigval.min();
    arma::vec beta = eigval - lambda_min;
    arma::uvec imin = arma::find(beta == 0.0);

    double C1 = 0.0, C2 = 0.0, C3 = 0.0;
    arma::uvec is_imin(gq.n_elem, arma::fill::zeros);
    for (arma::uword k = 0; k < imin.n_elem; k++) is_imin(imin(k)) = 1;
    for (arma::uword i = 0; i < gq.n_elem; i++) {
      C3 += gq(i) * gq(i);
      if (is_imin(i)) {
        C2 += gq(i) * gq(i);
      } else {
        double t = gq(i) / beta(i);
        C1 += t * t;
      }
    }

    if (C2 > 0.0 || C1 > r * r) {
      is_easy = true;
      is_hard = (C2 == 0.0);
      // ##### easy cases #####
      double root = trust_solve_lagrange_multiplier(beta, gq, r, C1, C2, C3);
      arma::vec wtry = gq / (beta + root);
      ptry = -(eigvec * wtry);
      steptype = is_hard ? TRUST_STEP_HARD_EASY : TRUST_STEP_EASY_EASY;
    } else {
      is_hard = true;
      is_easy = false;
      // ##### hard-hard case #####
      arma::vec wtry = gq / beta;  // NaN at imin (0/0); overwritten below,
                                    // exactly mirroring wtry[imin] <- 0 in R
      for (arma::uword k = 0; k < imin.n_elem; k++) wtry(imin(k)) = 0.0;
      ptry = -(eigvec * wtry);
      double utry_sq = r * r - arma::dot(ptry, ptry);
      if (utry_sq > 0.0) {
        double utry = std::sqrt(utry_sq);
        arma::vec vtry = eigvec.col(imin(0));
        ptry = ptry + utry * vtry;
      }
      steptype = TRUST_STEP_HARD_HARD;
    }
  }

  res.ptry = ptry;
  res.is_newton = is_newton;
  res.is_hard = is_hard;
  res.is_easy = is_easy;
  res.steptype = steptype;
  return res;
}
