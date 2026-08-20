#include "trust_rootfind.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <limits>

namespace {

// fred(beep) from upstream trust.R, verbatim translation including the
// beep == 0 special case (only ever invoked with beep == beta_dn == 0,
// i.e. the "hard-easy" boundary check).
double fred(double beep, const arma::vec &beta, const arma::vec &gq,
            double r, double C1, double C2) {
  if (beep == 0.0) {
    if (C2 > 0.0) return -1.0 / r;
    return std::sqrt(1.0 / C1) - 1.0 / r;
  }
  double s = 0.0;
  for (arma::uword i = 0; i < gq.n_elem; i++) {
    double w = gq(i) / (beta(i) + beep);
    s += w * w;
  }
  return std::sqrt(1.0 / s) - 1.0 / r;
}

}  // namespace

double trust_solve_lagrange_multiplier(const arma::vec &beta,
                                        const arma::vec &gq, double r,
                                        double C1, double C2, double C3) {
  double beta_dn = std::sqrt(C2) / r;
  double beta_up = std::sqrt(C3) / r;

  double f_up = fred(beta_up, beta, gq, r, C1, C2);
  if (f_up <= 0.0) return beta_up;

  double f_dn = fred(beta_dn, beta, gq, r, C1, C2);
  if (f_dn >= 0.0) return beta_dn;

  // f_dn < 0 < f_up: bracketed root, matching stats::uniroot()'s
  // requirement (and upstream's implicit assumption) exactly.
  auto f = [&](double beep) { return fred(beep, beta, gq, r, C1, C2); };

  // Match R's uniroot() default tol = .Machine$double.eps^0.25 on the
  // width of the bracketing interval.
  const double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
  auto tol_pred = [tol](double a, double b) { return std::fabs(b - a) <= tol; };

  boost::uintmax_t max_iter = 1000;
  std::pair<double, double> bracket = boost::math::tools::toms748_solve(
      f, beta_dn, beta_up, f_dn, f_up, tol_pred, max_iter);

  double fa = f(bracket.first);
  double fb = f(bracket.second);
  return std::fabs(fa) <= std::fabs(fb) ? bracket.first : bracket.second;
}
