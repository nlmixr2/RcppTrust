#ifndef TRUST_ROOTFIND_H
#define TRUST_ROOTFIND_H

#include <RcppArmadillo.h>

/* Port of the root-find used in the "easy" branch of upstream trust()'s
 * trust-region subproblem: solve
 *
 *   fred(beep) = sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r  ==  0
 *
 * for beep on [beta_dn, beta_up], where beta_dn = sqrt(C2) / r and
 * beta_up = sqrt(C3) / r. Upstream calls stats::uniroot(); this port
 * uses boost::math::tools::toms748_solve (both are derivative-free
 * bracketing methods), after applying the exact same boundary
 * short-circuits upstream applies before ever calling the root finder.
 *
 * beta: eigval - lambda_min, length d, all entries >= 0 (some may be 0).
 * gq: eigenvectors.t() * g, length d.
 * C1 = sum((gq / beta)[beta != 0]^2)
 * C2 = sum(gq[beta == 0]^2)
 * C3 = sum(gq^2)
 *
 * Returns the root (uout$root in the R code). */
double trust_solve_lagrange_multiplier(const arma::vec &beta,
                                        const arma::vec &gq, double r,
                                        double C1, double C2, double C3);

#endif /* TRUST_ROOTFIND_H */
