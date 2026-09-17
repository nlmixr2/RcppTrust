// Thread-safe C++17 port of M. J. D. Powell's NEWUOA, exactly as shipped
// (with R-specific edits by John Nash and Douglas Bates) in the CRAN
// package minqa 1.2.8: src/newuoa.f, newuob.f, trsapp.f, biglag.f,
// bigden.f and update.f.
//
// The translation is line by line: Powell's variable names, the partition
// of the workspace W, column-major 2-D arrays (accessed through 1-based
// index macros), statement labels (as goto targets) and the order of the
// floating point operations are all kept, so that results are bitwise
// identical to the Fortran compiled without FMA contraction.
//
// Differences from the Fortran, all of them required for thread safety:
//  * CALFUN(N,X,IPRINT) is the per-call MinqaCalfun object; MINQIT/MINQIR
//    are its member functions (no Rprintf, no global state).
//  * Local variables are initialized to zero (the Fortran leaves them
//    undefined, but never reads them before they are set).

#include <algorithm>
#include <cmath>
#include <vector>

#include "minqa_common.h"
#include "minqa_types.h"

// Forbid FMA contraction (a*b + c -> fma) so the port keeps the exact IEEE
// operation sequence of minqa's Fortran as built on x86_64 (and the
// bit-for-bit validation holds) on targets that contract by default, such
// as aarch64. Placed after the headers so only this file's code is affected.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

namespace {

// ---------------------------------------------------------------------------
// trsapp.f: SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
//                              D,G,HD,HS,CRVMIN)
// ---------------------------------------------------------------------------
void trsapp(int n, int npt, const double *xopt, const double *xpt,
            const double *gq, const double *hq, const double *pq,
            double delta, double *step, double *d, double *g, double *hd,
            double *hs, double &crvmin) {
#define XOPT(I) xopt[(I) - 1]
#define XPT(I, J) xpt[((I) - 1) + ((J) - 1) * npt]
#define GQ(I) gq[(I) - 1]
#define HQ(I) hq[(I) - 1]
#define PQ(I) pq[(I) - 1]
#define STEP(I) step[(I) - 1]
#define D(I) d[(I) - 1]
#define G(I) g[(I) - 1]
#define HD(I) hd[(I) - 1]
#define HS(I) hs[(I) - 1]
  double half, zero, twopi, delsq, qred, dd, ds, ss, gg, ggbeg, temp, bstep,
      dhd, alpha, qadd, ggsav, sg, shs, sgk, angtest, tempa, tempb, dg, dhs,
      cf, qbeg, qsav, qmin, angle, cth, sth, qnew, reduc, ratio;
  int iterc, itermax, itersw, i, j, k, ih, isave, iu;
  qred = dd = ds = ss = gg = ggbeg = temp = bstep = dhd = alpha = qadd =
      ggsav = sg = shs = sgk = angtest = tempa = tempb = dg = dhs = cf =
          qbeg = qsav = qmin = angle = cth = sth = qnew = reduc = ratio = 0.0;
  ih = isave = iu = 0;

  half = 0.5e0;
  zero = 0.0e0;
  twopi = 8.0e0 * std::atan(1.0e0);
  delsq = delta * delta;
  iterc = 0;
  itermax = n;
  itersw = itermax;
  for (i = 1; i <= n; i++) {
    D(i) = XOPT(i);
  }
  goto L170;

  // Prepare for the first line search.
L20:
  qred = zero;
  dd = zero;
  for (i = 1; i <= n; i++) {
    STEP(i) = zero;
    HS(i) = zero;
    G(i) = GQ(i) + HD(i);
    D(i) = -G(i);
    dd = dd + D(i) * D(i);
  }
  crvmin = zero;
  if (dd == zero) goto L160;
  ds = zero;
  ss = zero;
  gg = dd;
  ggbeg = gg;

  // Calculate the step to the trust region boundary and the product HD.
L40:
  iterc = iterc + 1;
  temp = delsq - ss;
  bstep = temp / (ds + std::sqrt(ds * ds + dd * temp));
  goto L170;
L50:
  dhd = zero;
  for (j = 1; j <= n; j++) {
    dhd = dhd + D(j) * HD(j);
  }

  // Update CRVMIN and set the step-length ALPHA.
  alpha = bstep;
  if (dhd > zero) {
    temp = dhd / dd;
    if (iterc == 1) crvmin = temp;
    crvmin = std::min(crvmin, temp);
    alpha = std::min(alpha, gg / dhd);
  }
  qadd = alpha * (gg - half * alpha * dhd);
  qred = qred + qadd;

  // Update STEP and HS.
  ggsav = gg;
  gg = zero;
  for (i = 1; i <= n; i++) {
    STEP(i) = STEP(i) + alpha * D(i);
    HS(i) = HS(i) + alpha * HD(i);
    gg = gg + (G(i) + HS(i)) * (G(i) + HS(i));
  }

  // Begin another conjugate direction iteration if required.
  if (alpha < bstep) {
    if (qadd <= 0.01e0 * qred) goto L160;
    if (gg <= 1.0e-4 * ggbeg) goto L160;
    if (iterc == itermax) goto L160;
    temp = gg / ggsav;
    dd = zero;
    ds = zero;
    ss = zero;
    for (i = 1; i <= n; i++) {
      D(i) = temp * D(i) - G(i) - HS(i);
      dd = dd + D(i) * D(i);
      ds = ds + D(i) * STEP(i);
      ss = ss + STEP(i) * STEP(i);
    }
    if (ds <= zero) goto L160;
    if (ss < delsq) goto L40;
  }
  crvmin = zero;
  itersw = iterc;

  // Test whether an alternative iteration is required.
L90:
  if (gg <= 1.0e-4 * ggbeg) goto L160;
  sg = zero;
  shs = zero;
  for (i = 1; i <= n; i++) {
    sg = sg + STEP(i) * G(i);
    shs = shs + STEP(i) * HS(i);
  }
  sgk = sg + shs;
  angtest = sgk / std::sqrt(gg * delsq);
  if (angtest <= -0.99e0) goto L160;

  // Begin the alternative iteration by calculating D and HD and some
  // scalar products.
  iterc = iterc + 1;
  temp = std::sqrt(delsq * gg - sgk * sgk);
  tempa = delsq / temp;
  tempb = sgk / temp;
  for (i = 1; i <= n; i++) {
    D(i) = tempa * (G(i) + HS(i)) - tempb * STEP(i);
  }
  goto L170;
L120:
  dg = zero;
  dhd = zero;
  dhs = zero;
  for (i = 1; i <= n; i++) {
    dg = dg + D(i) * G(i);
    dhd = dhd + HD(i) * D(i);
    dhs = dhs + HD(i) * STEP(i);
  }

  // Seek the value of the angle that minimizes Q.
  cf = half * (shs - dhd);
  qbeg = sg + cf;
  qsav = qbeg;
  qmin = qbeg;
  isave = 0;
  iu = 49;
  temp = twopi / static_cast<double>(iu + 1);
  for (i = 1; i <= iu; i++) {
    angle = static_cast<double>(i) * temp;
    cth = std::cos(angle);
    sth = std::sin(angle);
    qnew = (sg + cf * cth) * cth + (dg + dhs * cth) * sth;
    if (qnew < qmin) {
      qmin = qnew;
      isave = i;
      tempa = qsav;
    } else if (i == isave + 1) {
      tempb = qnew;
    }
    qsav = qnew;
  }
  if (static_cast<double>(isave) == zero) tempa = qnew;
  if (isave == iu) tempb = qbeg;
  angle = zero;
  if (tempa != tempb) {
    tempa = tempa - qmin;
    tempb = tempb - qmin;
    angle = half * (tempa - tempb) / (tempa + tempb);
  }
  angle = temp * (static_cast<double>(isave) + angle);

  // Calculate the new STEP and HS. Then test for convergence.
  cth = std::cos(angle);
  sth = std::sin(angle);
  reduc = qbeg - (sg + cf * cth) * cth - (dg + dhs * cth) * sth;
  gg = zero;
  for (i = 1; i <= n; i++) {
    STEP(i) = cth * STEP(i) + sth * D(i);
    HS(i) = cth * HS(i) + sth * HD(i);
    gg = gg + (G(i) + HS(i)) * (G(i) + HS(i));
  }
  qred = qred + reduc;
  ratio = reduc / qred;
  if (iterc < itermax && ratio > 0.01e0) goto L90;
L160:
  return;

  // The following instructions act as a subroutine for setting the vector
  // HD to the vector D multiplied by the second derivative matrix of Q.
  // They are called from three different places, which are distinguished
  // by the value of ITERC.
L170:
  for (i = 1; i <= n; i++) {
    HD(i) = zero;
  }
  for (k = 1; k <= npt; k++) {
    temp = zero;
    for (j = 1; j <= n; j++) {
      temp = temp + XPT(k, j) * D(j);
    }
    temp = temp * PQ(k);
    for (i = 1; i <= n; i++) {
      HD(i) = HD(i) + temp * XPT(k, i);
    }
  }
  ih = 0;
  for (j = 1; j <= n; j++) {
    for (i = 1; i <= j; i++) {
      ih = ih + 1;
      if (i < j) HD(j) = HD(j) + HQ(ih) * D(i);
      HD(i) = HD(i) + HQ(ih) * D(j);
    }
  }
  if (iterc == 0) goto L20;
  if (iterc <= itersw) goto L50;
  goto L120;
#undef XOPT
#undef XPT
#undef GQ
#undef HQ
#undef PQ
#undef STEP
#undef D
#undef G
#undef HD
#undef HS
}

// ---------------------------------------------------------------------------
// biglag.f: SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,
//                              DELTA,D,ALPHA,HCOL,GC,GD,S,W)
// ---------------------------------------------------------------------------
void biglag(int n, int npt, const double *xopt, const double *xpt,
            const double *bmat, const double *zmat, int idz, int ndim,
            int knew, double delta, double *d, double &alpha, double *hcol,
            double *gc, double *gd, double *s, double *w) {
#define XOPT(I) xopt[(I) - 1]
#define XPT(I, J) xpt[((I) - 1) + ((J) - 1) * npt]
#define BMAT(I, J) bmat[((I) - 1) + ((J) - 1) * ndim]
#define ZMAT(I, J) zmat[((I) - 1) + ((J) - 1) * npt]
#define D(I) d[(I) - 1]
#define HCOL(I) hcol[(I) - 1]
#define GC(I) gc[(I) - 1]
#define GD(I) gd[(I) - 1]
#define S(I) s[(I) - 1]
#define W(I) w[(I) - 1]
  double half, one, zero, twopi, delsq, temp, dd, sum, gg, sp, dhd, scale,
      tau, ss, denom, cf1, cf2, cf3, cf4, cf5, taubeg, taumax, tauold, angle,
      cth, sth, tempa, tempb, step;
  int nptm, iterc, i, j, k, isave, iu;
  temp = dd = sum = gg = sp = dhd = scale = tau = ss = denom = cf1 = cf2 =
      cf3 = cf4 = cf5 = taubeg = taumax = tauold = angle = cth = sth = tempa =
          tempb = step = 0.0;
  isave = iu = 0;

  half = 0.5e0;
  one = 1.0e0;
  zero = 0.0e0;
  twopi = 8.0e0 * std::atan(one);
  delsq = delta * delta;
  nptm = npt - n - 1;

  // Set the first NPT components of HCOL to the leading elements of the
  // KNEW-th column of H.
  iterc = 0;
  for (k = 1; k <= npt; k++) {
    HCOL(k) = zero;
  }
  for (j = 1; j <= nptm; j++) {
    temp = ZMAT(knew, j);
    if (j < idz) temp = -temp;
    for (k = 1; k <= npt; k++) {
      HCOL(k) = HCOL(k) + temp * ZMAT(k, j);
    }
  }
  alpha = HCOL(knew);

  // Set the unscaled initial direction D. Form the gradient of LFUNC at
  // XOPT, and multiply D by the second derivative matrix of LFUNC.
  dd = zero;
  for (i = 1; i <= n; i++) {
    D(i) = XPT(knew, i) - XOPT(i);
    GC(i) = BMAT(knew, i);
    GD(i) = zero;
    dd = dd + D(i) * D(i);
  }
  for (k = 1; k <= npt; k++) {
    temp = zero;
    sum = zero;
    for (j = 1; j <= n; j++) {
      temp = temp + XPT(k, j) * XOPT(j);
      sum = sum + XPT(k, j) * D(j);
    }
    temp = HCOL(k) * temp;
    sum = HCOL(k) * sum;
    for (i = 1; i <= n; i++) {
      GC(i) = GC(i) + temp * XPT(k, i);
      GD(i) = GD(i) + sum * XPT(k, i);
    }
  }

  // Scale D and GD, with a sign change if required. Set S to another
  // vector in the initial two dimensional subspace.
  gg = zero;
  sp = zero;
  dhd = zero;
  for (i = 1; i <= n; i++) {
    gg = gg + GC(i) * GC(i);
    sp = sp + D(i) * GC(i);
    dhd = dhd + D(i) * GD(i);
  }
  scale = delta / std::sqrt(dd);
  if (sp * dhd < zero) scale = -scale;
  temp = zero;
  if (sp * sp > 0.99e0 * dd * gg) temp = one;
  tau = scale * (std::fabs(sp) + half * scale * std::fabs(dhd));
  if (gg * delsq < 0.01e0 * tau * tau) temp = one;
  for (i = 1; i <= n; i++) {
    D(i) = scale * D(i);
    GD(i) = scale * GD(i);
    S(i) = GC(i) + temp * GD(i);
  }

  // Begin the iteration by overwriting S with a vector that has the
  // required length and direction, except that termination occurs if
  // the given D and S are nearly parallel.
L80:
  iterc = iterc + 1;
  dd = zero;
  sp = zero;
  ss = zero;
  for (i = 1; i <= n; i++) {
    dd = dd + D(i) * D(i);
    sp = sp + D(i) * S(i);
    ss = ss + S(i) * S(i);
  }
  temp = dd * ss - sp * sp;
  if (temp <= 1.0e-8 * dd * ss) goto L160;
  denom = std::sqrt(temp);
  for (i = 1; i <= n; i++) {
    S(i) = (dd * S(i) - sp * D(i)) / denom;
    W(i) = zero;
  }

  // Calculate the coefficients of the objective function on the circle,
  // beginning with the multiplication of S by the second derivative matrix.
  for (k = 1; k <= npt; k++) {
    sum = zero;
    for (j = 1; j <= n; j++) {
      sum = sum + XPT(k, j) * S(j);
    }
    sum = HCOL(k) * sum;
    for (i = 1; i <= n; i++) {
      W(i) = W(i) + sum * XPT(k, i);
    }
  }
  cf1 = zero;
  cf2 = zero;
  cf3 = zero;
  cf4 = zero;
  cf5 = zero;
  for (i = 1; i <= n; i++) {
    cf1 = cf1 + S(i) * W(i);
    cf2 = cf2 + D(i) * GC(i);
    cf3 = cf3 + S(i) * GC(i);
    cf4 = cf4 + D(i) * GD(i);
    cf5 = cf5 + S(i) * GD(i);
  }
  cf1 = half * cf1;
  cf4 = half * cf4 - cf1;

  // Seek the value of the angle that maximizes the modulus of TAU.
  taubeg = cf1 + cf2 + cf4;
  taumax = taubeg;
  tauold = taubeg;
  isave = 0;
  iu = 49;
  temp = twopi / static_cast<double>(iu + 1);
  for (i = 1; i <= iu; i++) {
    angle = static_cast<double>(i) * temp;
    cth = std::cos(angle);
    sth = std::sin(angle);
    tau = cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth) * sth;
    if (std::fabs(tau) > std::fabs(taumax)) {
      taumax = tau;
      isave = i;
      tempa = tauold;
    } else if (i == isave + 1) {
      tempb = tau;
    }
    tauold = tau;
  }
  if (isave == 0) tempa = tau;
  if (isave == iu) tempb = taubeg;
  step = zero;
  if (tempa != tempb) {
    tempa = tempa - taumax;
    tempb = tempb - taumax;
    step = half * (tempa - tempb) / (tempa + tempb);
  }
  angle = temp * (static_cast<double>(isave) + step);

  // Calculate the new D and GD. Then test for convergence.
  cth = std::cos(angle);
  sth = std::sin(angle);
  tau = cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth) * sth;
  for (i = 1; i <= n; i++) {
    D(i) = cth * D(i) + sth * S(i);
    GD(i) = cth * GD(i) + sth * W(i);
    S(i) = GC(i) + GD(i);
  }
  if (std::fabs(tau) <= 1.1e0 * std::fabs(taubeg)) goto L160;
  if (iterc < n) goto L80;
L160:
  return;
#undef XOPT
#undef XPT
#undef BMAT
#undef ZMAT
#undef D
#undef HCOL
#undef GC
#undef GD
#undef S
#undef W
}

// ---------------------------------------------------------------------------
// bigden.f: SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
//                              KNEW,D,W,VLAG,BETA,S,WVEC,PROD)
// ---------------------------------------------------------------------------
void bigden(int n, int npt, const double *xopt, const double *xpt,
            const double *bmat, const double *zmat, int idz, int ndim,
            int kopt, int knew, double *d, double *w, double *vlag,
            double &beta, double *s, double *wvec, double *prod) {
#define XOPT(I) xopt[(I) - 1]
#define XPT(I, J) xpt[((I) - 1) + ((J) - 1) * npt]
#define BMAT(I, J) bmat[((I) - 1) + ((J) - 1) * ndim]
#define ZMAT(I, J) zmat[((I) - 1) + ((J) - 1) * npt]
#define D(I) d[(I) - 1]
#define W(I) w[(I) - 1]
#define VLAG(I) vlag[(I) - 1]
#define S(I) s[(I) - 1]
#define WVEC(I, J) wvec[((I) - 1) + ((J) - 1) * ndim]
#define PROD(I, J) prod[((I) - 1) + ((J) - 1) * ndim]
#define DEN(I) den[(I) - 1]
#define DENEX(I) denex[(I) - 1]
#define PAR(I) par[(I) - 1]
  double den[9] = {0}, denex[9] = {0}, par[9] = {0};
  double half, one, quart, two, zero, twopi, temp, alpha, dd, ds, ss, xoptsq,
      dtest, dstemp, sstemp, diff, ssden, densav, xoptd, xopts, tempa, tempb,
      tempc, sum, denold, denmax, angle, sumold, step, tau;
  int nptm, i, j, k, ksav, iterc, ip, jc, nw, isave, iu;
  temp = alpha = dd = ds = ss = xoptsq = dtest = dstemp = sstemp = diff =
      ssden = densav = xoptd = xopts = tempa = tempb = tempc = sum = denold =
          denmax = angle = sumold = step = tau = 0.0;
  ksav = iterc = ip = jc = nw = isave = iu = 0;

  half = 0.5e0;
  one = 1.0e0;
  quart = 0.25e0;
  two = 2.0e0;
  zero = 0.0e0;
  twopi = 8.0e0 * std::atan(one);
  nptm = npt - n - 1;

  // Store the first NPT elements of the KNEW-th column of H in W(N+1)
  // to W(N+NPT).
  for (k = 1; k <= npt; k++) {
    W(n + k) = zero;
  }
  for (j = 1; j <= nptm; j++) {
    temp = ZMAT(knew, j);
    if (j < idz) temp = -temp;
    for (k = 1; k <= npt; k++) {
      W(n + k) = W(n + k) + temp * ZMAT(k, j);
    }
  }
  alpha = W(n + knew);

  // The initial search direction D is taken from the last call of BIGLAG,
  // and the initial S is set below, usually to the direction from X_OPT
  // to X_KNEW, but a different direction to an interpolation point may
  // be chosen, in order to prevent S from being nearly parallel to D.
  dd = zero;
  ds = zero;
  ss = zero;
  xoptsq = zero;
  for (i = 1; i <= n; i++) {
    dd = dd + D(i) * D(i);
    S(i) = XPT(knew, i) - XOPT(i);
    ds = ds + D(i) * S(i);
    ss = ss + S(i) * S(i);
    xoptsq = xoptsq + XOPT(i) * XOPT(i);
  }
  if (ds * ds > 0.99e0 * dd * ss) {
    ksav = knew;
    dtest = ds * ds / ss;
    for (k = 1; k <= npt; k++) {
      if (k != kopt) {
        dstemp = zero;
        sstemp = zero;
        for (i = 1; i <= n; i++) {
          diff = XPT(k, i) - XOPT(i);
          dstemp = dstemp + D(i) * diff;
          sstemp = sstemp + diff * diff;
        }
        if (dstemp * dstemp / sstemp < dtest) {
          ksav = k;
          dtest = dstemp * dstemp / sstemp;
          ds = dstemp;
          ss = sstemp;
        }
      }
    }
    for (i = 1; i <= n; i++) {
      S(i) = XPT(ksav, i) - XOPT(i);
    }
  }
  ssden = dd * ss - ds * ds;
  iterc = 0;
  densav = zero;

  // Begin the iteration by overwriting S with a vector that has the
  // required length and direction.
L70:
  iterc = iterc + 1;
  temp = one / std::sqrt(ssden);
  xoptd = zero;
  xopts = zero;
  for (i = 1; i <= n; i++) {
    S(i) = temp * (dd * S(i) - ds * D(i));
    xoptd = xoptd + XOPT(i) * D(i);
    xopts = xopts + XOPT(i) * S(i);
  }

  // Set the coefficients of the first two terms of BETA.
  tempa = half * xoptd * xoptd;
  tempb = half * xopts * xopts;
  DEN(1) = dd * (xoptsq + half * dd) + tempa + tempb;
  DEN(2) = two * xoptd * dd;
  DEN(3) = two * xopts * dd;
  DEN(4) = tempa - tempb;
  DEN(5) = xoptd * xopts;
  for (i = 6; i <= 9; i++) {
    DEN(i) = zero;
  }

  // Put the coefficients of Wcheck in WVEC.
  for (k = 1; k <= npt; k++) {
    tempa = zero;
    tempb = zero;
    tempc = zero;
    for (i = 1; i <= n; i++) {
      tempa = tempa + XPT(k, i) * D(i);
      tempb = tempb + XPT(k, i) * S(i);
      tempc = tempc + XPT(k, i) * XOPT(i);
    }
    WVEC(k, 1) = quart * (tempa * tempa + tempb * tempb);
    WVEC(k, 2) = tempa * tempc;
    WVEC(k, 3) = tempb * tempc;
    WVEC(k, 4) = quart * (tempa * tempa - tempb * tempb);
    WVEC(k, 5) = half * tempa * tempb;
  }
  for (i = 1; i <= n; i++) {
    ip = i + npt;
    WVEC(ip, 1) = zero;
    WVEC(ip, 2) = D(i);
    WVEC(ip, 3) = S(i);
    WVEC(ip, 4) = zero;
    WVEC(ip, 5) = zero;
  }

  // Put the coefficents of THETA*Wcheck in PROD.
  for (jc = 1; jc <= 5; jc++) {
    nw = npt;
    if (jc == 2 || jc == 3) nw = ndim;
    for (k = 1; k <= npt; k++) {
      PROD(k, jc) = zero;
    }
    for (j = 1; j <= nptm; j++) {
      sum = zero;
      for (k = 1; k <= npt; k++) {
        sum = sum + ZMAT(k, j) * WVEC(k, jc);
      }
      if (j < idz) sum = -sum;
      for (k = 1; k <= npt; k++) {
        PROD(k, jc) = PROD(k, jc) + sum * ZMAT(k, j);
      }
    }
    if (nw == ndim) {
      for (k = 1; k <= npt; k++) {
        sum = zero;
        for (j = 1; j <= n; j++) {
          sum = sum + BMAT(k, j) * WVEC(npt + j, jc);
        }
        PROD(k, jc) = PROD(k, jc) + sum;
      }
    }
    for (j = 1; j <= n; j++) {
      sum = zero;
      for (i = 1; i <= nw; i++) {
        sum = sum + BMAT(i, j) * WVEC(i, jc);
      }
      PROD(npt + j, jc) = sum;
    }
  }

  // Include in DEN the part of BETA that depends on THETA.
  for (k = 1; k <= ndim; k++) {
    sum = zero;
    for (i = 1; i <= 5; i++) {
      PAR(i) = half * PROD(k, i) * WVEC(k, i);
      sum = sum + PAR(i);
    }
    DEN(1) = DEN(1) - PAR(1) - sum;
    tempa = PROD(k, 1) * WVEC(k, 2) + PROD(k, 2) * WVEC(k, 1);
    tempb = PROD(k, 2) * WVEC(k, 4) + PROD(k, 4) * WVEC(k, 2);
    tempc = PROD(k, 3) * WVEC(k, 5) + PROD(k, 5) * WVEC(k, 3);
    DEN(2) = DEN(2) - tempa - half * (tempb + tempc);
    DEN(6) = DEN(6) - half * (tempb - tempc);
    tempa = PROD(k, 1) * WVEC(k, 3) + PROD(k, 3) * WVEC(k, 1);
    tempb = PROD(k, 2) * WVEC(k, 5) + PROD(k, 5) * WVEC(k, 2);
    tempc = PROD(k, 3) * WVEC(k, 4) + PROD(k, 4) * WVEC(k, 3);
    DEN(3) = DEN(3) - tempa - half * (tempb - tempc);
    DEN(7) = DEN(7) - half * (tempb + tempc);
    tempa = PROD(k, 1) * WVEC(k, 4) + PROD(k, 4) * WVEC(k, 1);
    DEN(4) = DEN(4) - tempa - PAR(2) + PAR(3);
    tempa = PROD(k, 1) * WVEC(k, 5) + PROD(k, 5) * WVEC(k, 1);
    tempb = PROD(k, 2) * WVEC(k, 3) + PROD(k, 3) * WVEC(k, 2);
    DEN(5) = DEN(5) - tempa - half * tempb;
    DEN(8) = DEN(8) - PAR(4) + PAR(5);
    tempa = PROD(k, 4) * WVEC(k, 5) + PROD(k, 5) * WVEC(k, 4);
    DEN(9) = DEN(9) - half * tempa;
  }

  // Extend DEN so that it holds all the coefficients of DENOM.
  sum = zero;
  for (i = 1; i <= 5; i++) {
    PAR(i) = half * PROD(knew, i) * PROD(knew, i);
    sum = sum + PAR(i);
  }
  DENEX(1) = alpha * DEN(1) + PAR(1) + sum;
  tempa = two * PROD(knew, 1) * PROD(knew, 2);
  tempb = PROD(knew, 2) * PROD(knew, 4);
  tempc = PROD(knew, 3) * PROD(knew, 5);
  DENEX(2) = alpha * DEN(2) + tempa + tempb + tempc;
  DENEX(6) = alpha * DEN(6) + tempb - tempc;
  tempa = two * PROD(knew, 1) * PROD(knew, 3);
  tempb = PROD(knew, 2) * PROD(knew, 5);
  tempc = PROD(knew, 3) * PROD(knew, 4);
  DENEX(3) = alpha * DEN(3) + tempa + tempb - tempc;
  DENEX(7) = alpha * DEN(7) + tempb + tempc;
  tempa = two * PROD(knew, 1) * PROD(knew, 4);
  DENEX(4) = alpha * DEN(4) + tempa + PAR(2) - PAR(3);
  tempa = two * PROD(knew, 1) * PROD(knew, 5);
  DENEX(5) = alpha * DEN(5) + tempa + PROD(knew, 2) * PROD(knew, 3);
  DENEX(8) = alpha * DEN(8) + PAR(4) - PAR(5);
  DENEX(9) = alpha * DEN(9) + PROD(knew, 4) * PROD(knew, 5);

  // Seek the value of the angle that maximizes the modulus of DENOM.
  sum = DENEX(1) + DENEX(2) + DENEX(4) + DENEX(6) + DENEX(8);
  denold = sum;
  denmax = sum;
  isave = 0;
  iu = 49;
  temp = twopi / static_cast<double>(iu + 1);
  PAR(1) = one;
  for (i = 1; i <= iu; i++) {
    angle = static_cast<double>(i) * temp;
    PAR(2) = std::cos(angle);
    PAR(3) = std::sin(angle);
    for (j = 4; j <= 8; j += 2) {
      PAR(j) = PAR(2) * PAR(j - 2) - PAR(3) * PAR(j - 1);
      PAR(j + 1) = PAR(2) * PAR(j - 1) + PAR(3) * PAR(j - 2);
    }
    sumold = sum;
    sum = zero;
    for (j = 1; j <= 9; j++) {
      sum = sum + DENEX(j) * PAR(j);
    }
    if (std::fabs(sum) > std::fabs(denmax)) {
      denmax = sum;
      isave = i;
      tempa = sumold;
    } else if (i == isave + 1) {
      tempb = sum;
    }
  }
  if (isave == 0) tempa = sum;
  if (isave == iu) tempb = denold;
  step = zero;
  if (tempa != tempb) {
    tempa = tempa - denmax;
    tempb = tempb - denmax;
    step = half * (tempa - tempb) / (tempa + tempb);
  }
  angle = temp * (static_cast<double>(isave) + step);

  // Calculate the new parameters of the denominator, the new VLAG vector
  // and the new D. Then test for convergence.
  PAR(2) = std::cos(angle);
  PAR(3) = std::sin(angle);
  for (j = 4; j <= 8; j += 2) {
    PAR(j) = PAR(2) * PAR(j - 2) - PAR(3) * PAR(j - 1);
    PAR(j + 1) = PAR(2) * PAR(j - 1) + PAR(3) * PAR(j - 2);
  }
  beta = zero;
  denmax = zero;
  for (j = 1; j <= 9; j++) {
    beta = beta + DEN(j) * PAR(j);
    denmax = denmax + DENEX(j) * PAR(j);
  }
  for (k = 1; k <= ndim; k++) {
    VLAG(k) = zero;
    for (j = 1; j <= 5; j++) {
      VLAG(k) = VLAG(k) + PROD(k, j) * PAR(j);
    }
  }
  tau = VLAG(knew);
  dd = zero;
  tempa = zero;
  tempb = zero;
  for (i = 1; i <= n; i++) {
    D(i) = PAR(2) * D(i) + PAR(3) * S(i);
    W(i) = XOPT(i) + D(i);
    dd = dd + D(i) * D(i);
    tempa = tempa + D(i) * W(i);
    tempb = tempb + W(i) * W(i);
  }
  if (iterc >= n) goto L340;
  if (iterc > 1) densav = std::max(densav, denold);
  if (std::fabs(denmax) <= 1.1e0 * std::fabs(densav)) goto L340;
  densav = denmax;

  // Set S to half the gradient of the denominator with respect to D.
  // Then branch for the next iteration.
  for (i = 1; i <= n; i++) {
    temp = tempa * XOPT(i) + tempb * D(i) - VLAG(npt + i);
    S(i) = tau * BMAT(knew, i) + alpha * temp;
  }
  for (k = 1; k <= npt; k++) {
    sum = zero;
    for (j = 1; j <= n; j++) {
      sum = sum + XPT(k, j) * W(j);
    }
    temp = (tau * W(n + k) - alpha * VLAG(k)) * sum;
    for (i = 1; i <= n; i++) {
      S(i) = S(i) + temp * XPT(k, i);
    }
  }
  ss = zero;
  ds = zero;
  for (i = 1; i <= n; i++) {
    ss = ss + S(i) * S(i);
    ds = ds + D(i) * S(i);
  }
  ssden = dd * ss - ds * ds;
  if (ssden >= 1.0e-8 * dd * ss) goto L70;

  // Set the vector W before the RETURN from the subroutine.
L340:
  for (k = 1; k <= ndim; k++) {
    W(k) = zero;
    for (j = 1; j <= 5; j++) {
      W(k) = W(k) + WVEC(k, j) * PAR(j);
    }
  }
  VLAG(kopt) = VLAG(kopt) + one;
  return;
#undef XOPT
#undef XPT
#undef BMAT
#undef ZMAT
#undef D
#undef W
#undef VLAG
#undef S
#undef WVEC
#undef PROD
#undef DEN
#undef DENEX
#undef PAR
}

// ---------------------------------------------------------------------------
// update.f: SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
// ---------------------------------------------------------------------------
void update(int n, int npt, double *bmat, double *zmat, int &idz, int ndim,
            double *vlag, double beta, int knew, double *w) {
#define BMAT(I, J) bmat[((I) - 1) + ((J) - 1) * ndim]
#define ZMAT(I, J) zmat[((I) - 1) + ((J) - 1) * npt]
#define VLAG(I) vlag[(I) - 1]
#define W(I) w[(I) - 1]
  double one, zero, temp, tempa, tempb, alpha, tau, tausq, denom, scala,
      scalb;
  int nptm, jl, i, j, iflag, ja, jb, jp;
  temp = tempa = tempb = alpha = tau = tausq = denom = scala = scalb = 0.0;
  ja = jb = jp = 0;

  one = 1.0e0;
  zero = 0.0e0;
  nptm = npt - n - 1;

  // Apply the rotations that put zeros in the KNEW-th row of ZMAT.
  jl = 1;
  for (j = 2; j <= nptm; j++) {
    if (j == idz) {
      jl = idz;
    } else if (ZMAT(knew, j) != zero) {
      temp = std::sqrt(ZMAT(knew, jl) * ZMAT(knew, jl) +
                       ZMAT(knew, j) * ZMAT(knew, j));
      tempa = ZMAT(knew, jl) / temp;
      tempb = ZMAT(knew, j) / temp;
      for (i = 1; i <= npt; i++) {
        temp = tempa * ZMAT(i, jl) + tempb * ZMAT(i, j);
        ZMAT(i, j) = tempa * ZMAT(i, j) - tempb * ZMAT(i, jl);
        ZMAT(i, jl) = temp;
      }
      ZMAT(knew, j) = zero;
    }
  }

  // Put the first NPT components of the KNEW-th column of HLAG into W,
  // and calculate the parameters of the updating formula.
  tempa = ZMAT(knew, 1);
  if (idz >= 2) tempa = -tempa;
  if (jl > 1) tempb = ZMAT(knew, jl);
  for (i = 1; i <= npt; i++) {
    W(i) = tempa * ZMAT(i, 1);
    if (jl > 1) W(i) = W(i) + tempb * ZMAT(i, jl);
  }
  alpha = W(knew);
  tau = VLAG(knew);
  tausq = tau * tau;
  denom = alpha * beta + tausq;
  VLAG(knew) = VLAG(knew) - one;

  // Complete the updating of ZMAT when there is only one nonzero element
  // in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
  // then the first column of ZMAT will be exchanged with another one later.
  iflag = 0;
  if (jl == 1) {
    temp = std::sqrt(std::fabs(denom));
    tempb = tempa / temp;
    tempa = tau / temp;
    for (i = 1; i <= npt; i++) {
      ZMAT(i, 1) = tempa * ZMAT(i, 1) - tempb * VLAG(i);
    }
    if (idz == 1 && temp < zero) idz = 2;
    if (idz >= 2 && temp >= zero) iflag = 1;
  } else {
    // Complete the updating of ZMAT in the alternative case.
    ja = 1;
    if (beta >= zero) ja = jl;
    jb = jl + 1 - ja;
    temp = ZMAT(knew, jb) / denom;
    tempa = temp * beta;
    tempb = temp * tau;
    temp = ZMAT(knew, ja);
    scala = one / std::sqrt(std::fabs(beta) * temp * temp + tausq);
    scalb = scala * std::sqrt(std::fabs(denom));
    for (i = 1; i <= npt; i++) {
      ZMAT(i, ja) = scala * (tau * ZMAT(i, ja) - temp * VLAG(i));
      ZMAT(i, jb) = scalb * (ZMAT(i, jb) - tempa * W(i) - tempb * VLAG(i));
    }
    if (denom <= zero) {
      if (beta < zero) idz = idz + 1;
      if (beta >= zero) iflag = 1;
    }
  }

  // IDZ is reduced in the following case, and usually the first column
  // of ZMAT is exchanged with a later one.
  if (iflag == 1) {
    idz = idz - 1;
    for (i = 1; i <= npt; i++) {
      temp = ZMAT(i, 1);
      ZMAT(i, 1) = ZMAT(i, idz);
      ZMAT(i, idz) = temp;
    }
  }

  // Finally, update the matrix BMAT.
  for (j = 1; j <= n; j++) {
    jp = npt + j;
    W(jp) = BMAT(knew, j);
    tempa = (alpha * VLAG(jp) - tau * W(jp)) / denom;
    tempb = (-beta * W(jp) - tau * VLAG(jp)) / denom;
    for (i = 1; i <= jp; i++) {
      BMAT(i, j) = BMAT(i, j) + tempa * VLAG(i) + tempb * W(i);
      if (i > npt) BMAT(jp, i - npt) = BMAT(i, j);
    }
  }
  return;
#undef BMAT
#undef ZMAT
#undef VLAG
#undef W
}

// ---------------------------------------------------------------------------
// newuob.f: SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
//             XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,IERR)
// ---------------------------------------------------------------------------
void newuob(int n, int npt, double *x, double rhobeg, double rhoend,
            int iprint, int maxfun, double *xbase, double *xopt,
            double *xnew, double *xpt, double *fval, double *gq, double *hq,
            double *pq, double *bmat, double *zmat, int ndim, double *d,
            double *vlag, double *w, int &ierr, MinqaCalfun &calfun) {
#define X(I) x[(I) - 1]
#define XBASE(I) xbase[(I) - 1]
#define XOPT(I) xopt[(I) - 1]
#define XNEW(I) xnew[(I) - 1]
#define XPT(I, J) xpt[((I) - 1) + ((J) - 1) * npt]
#define FVAL(I) fval[(I) - 1]
#define GQ(I) gq[(I) - 1]
#define HQ(I) hq[(I) - 1]
#define PQ(I) pq[(I) - 1]
#define BMAT(I, J) bmat[((I) - 1) + ((J) - 1) * ndim]
#define ZMAT(I, J) zmat[((I) - 1) + ((J) - 1) * npt]
#define D(I) d[(I) - 1]
#define VLAG(I) vlag[(I) - 1]
#define W(I) w[(I) - 1]
  double half, one, tenth, zero, rhosq, recip, reciq, xipt, xjpt, f, fbeg,
      fopt, temp, rho, delta, diffa, diffb, diffc, xoptsq, crvmin, dsq, dnorm,
      ratio, tempq, sum, sumz, dstep, alpha, suma, sumb, beta, bsum, dx,
      vquad, diff, fsave, detrat, hdiag, distsq, gqsq, gisq;
  int np, nh, nptm, nftest, i, j, k, ih, nf, nfm, nfmm, itemp, jpt, ipt,
      kopt, idz, itest, nfsav, knew, ip, jp, ksave, ktemp;
  xipt = xjpt = f = fbeg = fopt = temp = rho = delta = diffa = diffb = diffc =
      xoptsq = crvmin = dsq = dnorm = ratio = tempq = sum = sumz = dstep =
          alpha = suma = sumb = beta = bsum = dx = vquad = diff = fsave =
              detrat = hdiag = distsq = gqsq = gisq = 0.0;
  itemp = jpt = ipt = kopt = idz = itest = nfsav = knew = ip = jp = ksave =
      ktemp = 0;

  // Set some constants.
  half = 0.5e0;
  one = 1.0e0;
  tenth = 0.1e0;
  zero = 0.0e0;
  np = n + 1;
  nh = (n * np) / 2;
  nptm = npt - np;
  nftest = std::max(maxfun, 1);

  // Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
  for (j = 1; j <= n; j++) {
    XBASE(j) = X(j);
    for (k = 1; k <= npt; k++) {
      XPT(k, j) = zero;
    }
    for (i = 1; i <= ndim; i++) {
      BMAT(i, j) = zero;
    }
  }
  for (ih = 1; ih <= nh; ih++) {
    HQ(ih) = zero;
  }
  for (k = 1; k <= npt; k++) {
    PQ(k) = zero;
    for (j = 1; j <= nptm; j++) {
      ZMAT(k, j) = zero;
    }
  }

  // Begin the initialization procedure. NF becomes one more than the number
  // of function values so far. The coordinates of the displacement of the
  // next initial interpolation point from XBASE are set in XPT(NF,.).
  rhosq = rhobeg * rhobeg;
  recip = one / rhosq;
  reciq = std::sqrt(half) / rhosq;
  nf = 0;
L50:
  nfm = nf;
  nfmm = nf - n;
  nf = nf + 1;
  if (nfm <= 2 * n) {
    if (nfm >= 1 && nfm <= n) {
      XPT(nf, nfm) = rhobeg;
    } else if (nfm > n) {
      XPT(nf, nfmm) = -rhobeg;
    }
  } else {
    itemp = (nfmm - 1) / n;
    jpt = nfm - itemp * n - n;
    ipt = jpt + itemp;
    if (ipt > n) {
      itemp = jpt;
      jpt = ipt - n;
      ipt = itemp;
    }
    xipt = rhobeg;
    if (FVAL(ipt + np) < FVAL(ipt + 1)) xipt = -xipt;
    xjpt = rhobeg;
    if (FVAL(jpt + np) < FVAL(jpt + 1)) xjpt = -xjpt;
    XPT(nf, ipt) = xipt;
    XPT(nf, jpt) = xjpt;
  }

  // Calculate the next value of F, label 70 being reached immediately
  // after this calculation. The least function value so far and its index
  // are required.
  for (j = 1; j <= n; j++) {
    X(j) = XPT(nf, j) + XBASE(j);
  }
  goto L310;
L70:
  FVAL(nf) = f;
  if (nf == 1) {
    fbeg = f;
    fopt = f;
    kopt = 1;
  } else if (f < fopt) {
    fopt = f;
    kopt = nf;
  }

  // Set the nonzero initial elements of BMAT and the quadratic model in
  // the cases when NF is at most 2*N+1.
  if (nfm <= 2 * n) {
    if (nfm >= 1 && nfm <= n) {
      GQ(nfm) = (f - fbeg) / rhobeg;
      if (npt < nf + n) {
        BMAT(1, nfm) = -one / rhobeg;
        BMAT(nf, nfm) = one / rhobeg;
        BMAT(npt + nfm, nfm) = -half * rhosq;
      }
    } else if (nfm > n) {
      BMAT(nf - n, nfmm) = half / rhobeg;
      BMAT(nf, nfmm) = -half / rhobeg;
      ZMAT(1, nfmm) = -reciq - reciq;
      ZMAT(nf - n, nfmm) = reciq;
      ZMAT(nf, nfmm) = reciq;
      ih = (nfmm * (nfmm + 1)) / 2;
      temp = (fbeg - f) / rhobeg;
      HQ(ih) = (GQ(nfmm) - temp) / rhobeg;
      GQ(nfmm) = half * (GQ(nfmm) + temp);
    }

    // Set the off-diagonal second derivatives of the Lagrange functions and
    // the initial quadratic model.
  } else {
    ih = (ipt * (ipt - 1)) / 2 + jpt;
    if (xipt < zero) ipt = ipt + n;
    if (xjpt < zero) jpt = jpt + n;
    ZMAT(1, nfmm) = recip;
    ZMAT(nf, nfmm) = recip;
    ZMAT(ipt + 1, nfmm) = -recip;
    ZMAT(jpt + 1, nfmm) = -recip;
    HQ(ih) = (fbeg - FVAL(ipt + 1) - FVAL(jpt + 1) + f) / (xipt * xjpt);
  }
  if (nf < npt) goto L50;

  // Begin the iterative procedure, because the initial model is complete.
  rho = rhobeg;
  delta = rho;
  idz = 1;
  diffa = zero;
  diffb = zero;
  itest = 0;
  xoptsq = zero;
  for (i = 1; i <= n; i++) {
    XOPT(i) = XPT(kopt, i);
    xoptsq = xoptsq + XOPT(i) * XOPT(i);
  }
L90:
  nfsav = nf;

  // Generate the next trust region step and test its length. Set KNEW
  // to -1 if the purpose of the next F will be to improve the model.
L100:
  knew = 0;
  trsapp(n, npt, xopt, xpt, gq, hq, pq, delta, d, w, &W(np), &W(np + n),
         &W(np + 2 * n), crvmin);
  dsq = zero;
  for (i = 1; i <= n; i++) {
    dsq = dsq + D(i) * D(i);
  }
  dnorm = std::min(delta, std::sqrt(dsq));
  if (dnorm < half * rho) {
    knew = -1;
    delta = tenth * delta;
    ratio = -1.0e0;
    if (delta <= 1.5e0 * rho) delta = rho;
    if (nf <= nfsav + 2) goto L460;
    temp = 0.125e0 * crvmin * rho * rho;
    if (temp <= std::max(std::max(diffa, diffb), diffc)) goto L460;
    goto L490;
  }

  // Shift XBASE if XOPT may be too far from XBASE. First make the changes
  // to BMAT that do not depend on ZMAT.
L120:
  if (dsq <= 1.0e-3 * xoptsq) {
    tempq = 0.25e0 * xoptsq;
    for (k = 1; k <= npt; k++) {
      sum = zero;
      for (i = 1; i <= n; i++) {
        sum = sum + XPT(k, i) * XOPT(i);
      }
      temp = PQ(k) * sum;
      sum = sum - half * xoptsq;
      W(npt + k) = sum;
      for (i = 1; i <= n; i++) {
        GQ(i) = GQ(i) + temp * XPT(k, i);
        XPT(k, i) = XPT(k, i) - half * XOPT(i);
        VLAG(i) = BMAT(k, i);
        W(i) = sum * XPT(k, i) + tempq * XOPT(i);
        ip = npt + i;
        for (j = 1; j <= i; j++) {
          BMAT(ip, j) = BMAT(ip, j) + VLAG(i) * W(j) + W(i) * VLAG(j);
        }
      }
    }

    // Then the revisions of BMAT that depend on ZMAT are calculated.
    for (k = 1; k <= nptm; k++) {
      sumz = zero;
      for (i = 1; i <= npt; i++) {
        sumz = sumz + ZMAT(i, k);
        W(i) = W(npt + i) * ZMAT(i, k);
      }
      for (j = 1; j <= n; j++) {
        sum = tempq * sumz * XOPT(j);
        for (i = 1; i <= npt; i++) {
          sum = sum + W(i) * XPT(i, j);
        }
        VLAG(j) = sum;
        if (k < idz) sum = -sum;
        for (i = 1; i <= npt; i++) {
          BMAT(i, j) = BMAT(i, j) + sum * ZMAT(i, k);
        }
      }
      for (i = 1; i <= n; i++) {
        ip = i + npt;
        temp = VLAG(i);
        if (k < idz) temp = -temp;
        for (j = 1; j <= i; j++) {
          BMAT(ip, j) = BMAT(ip, j) + temp * VLAG(j);
        }
      }
    }

    // The following instructions complete the shift of XBASE, including
    // the changes to the parameters of the quadratic model.
    ih = 0;
    for (j = 1; j <= n; j++) {
      W(j) = zero;
      for (k = 1; k <= npt; k++) {
        W(j) = W(j) + PQ(k) * XPT(k, j);
        XPT(k, j) = XPT(k, j) - half * XOPT(j);
      }
      for (i = 1; i <= j; i++) {
        ih = ih + 1;
        if (i < j) GQ(j) = GQ(j) + HQ(ih) * XOPT(i);
        GQ(i) = GQ(i) + HQ(ih) * XOPT(j);
        HQ(ih) = HQ(ih) + W(i) * XOPT(j) + XOPT(i) * W(j);
        BMAT(npt + i, j) = BMAT(npt + j, i);
      }
    }
    for (j = 1; j <= n; j++) {
      XBASE(j) = XBASE(j) + XOPT(j);
      XOPT(j) = zero;
    }
    xoptsq = zero;
  }

  // Pick the model step if KNEW is positive. A different choice of D
  // may be made later, if the choice of D by BIGLAG causes substantial
  // cancellation in DENOM.
  if (knew > 0) {
    biglag(n, npt, xopt, xpt, bmat, zmat, idz, ndim, knew, dstep, d, alpha,
           vlag, &VLAG(npt + 1), w, &W(np), &W(np + n));
  }

  // Calculate VLAG and BETA for the current choice of D. The first NPT
  // components of W_check will be held in W.
  for (k = 1; k <= npt; k++) {
    suma = zero;
    sumb = zero;
    sum = zero;
    for (j = 1; j <= n; j++) {
      suma = suma + XPT(k, j) * D(j);
      sumb = sumb + XPT(k, j) * XOPT(j);
      sum = sum + BMAT(k, j) * D(j);
    }
    W(k) = suma * (half * suma + sumb);
    VLAG(k) = sum;
  }
  beta = zero;
  for (k = 1; k <= nptm; k++) {
    sum = zero;
    for (i = 1; i <= npt; i++) {
      sum = sum + ZMAT(i, k) * W(i);
    }
    if (k < idz) {
      beta = beta + sum * sum;
      sum = -sum;
    } else {
      beta = beta - sum * sum;
    }
    for (i = 1; i <= npt; i++) {
      VLAG(i) = VLAG(i) + sum * ZMAT(i, k);
    }
  }
  bsum = zero;
  dx = zero;
  for (j = 1; j <= n; j++) {
    sum = zero;
    for (i = 1; i <= npt; i++) {
      sum = sum + W(i) * BMAT(i, j);
    }
    bsum = bsum + sum * D(j);
    jp = npt + j;
    for (k = 1; k <= n; k++) {
      sum = sum + BMAT(jp, k) * D(k);
    }
    VLAG(jp) = sum;
    bsum = bsum + sum * D(j);
    dx = dx + D(j) * XOPT(j);
  }
  beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
  VLAG(kopt) = VLAG(kopt) + one;

  // If KNEW is positive and if the cancellation in DENOM is unacceptable,
  // then BIGDEN calculates an alternative model step, XNEW being used for
  // working space.
  if (knew > 0) {
    temp = one + alpha * beta / (VLAG(knew) * VLAG(knew));
    if (std::fabs(temp) <= 0.8e0) {
      bigden(n, npt, xopt, xpt, bmat, zmat, idz, ndim, kopt, knew, d, w, vlag,
             beta, xnew, &W(ndim + 1), &W(6 * ndim + 1));
    }
  }

  // Calculate the next value of the objective function.
L290:
  for (i = 1; i <= n; i++) {
    XNEW(i) = XOPT(i) + D(i);
    X(i) = XBASE(i) + XNEW(i);
  }
  nf = nf + 1;
L310:
  if (nf > nftest) {
    nf = nf - 1;
    ierr = 390;
    goto L530;
  }
  f = calfun(x);
  if (nf <= npt) goto L70;
  if (knew == -1) goto L530;

  // Use the quadratic model to predict the change in F due to the step D,
  // and set DIFF to the error of this prediction.
  vquad = zero;
  ih = 0;
  for (j = 1; j <= n; j++) {
    vquad = vquad + D(j) * GQ(j);
    for (i = 1; i <= j; i++) {
      ih = ih + 1;
      temp = D(i) * XNEW(j) + D(j) * XOPT(i);
      if (i == j) temp = half * temp;
      vquad = vquad + temp * HQ(ih);
    }
  }
  for (k = 1; k <= npt; k++) {
    vquad = vquad + PQ(k) * W(k);
  }
  diff = f - fopt - vquad;
  diffc = diffb;
  diffb = diffa;
  diffa = std::fabs(diff);
  if (dnorm > rho) nfsav = nf;

  // Update FOPT and XOPT if the new F is the least value of the objective
  // function so far. The branch when KNEW is positive occurs if D is not
  // a trust region step.
  fsave = fopt;
  if (f < fopt) {
    fopt = f;
    xoptsq = zero;
    for (i = 1; i <= n; i++) {
      XOPT(i) = XNEW(i);
      xoptsq = xoptsq + XOPT(i) * XOPT(i);
    }
  }
  ksave = knew;
  if (knew > 0) goto L410;

  // Pick the next value of DELTA after a trust region step.
  if (vquad >= zero) {
    ierr = 3701;
    goto L530;
  }
  ratio = (f - fsave) / vquad;
  if (ratio <= tenth) {
    delta = half * dnorm;
  } else if (ratio <= 0.7e0) {
    delta = std::max(half * delta, dnorm);
  } else {
    delta = std::max(half * delta, dnorm + dnorm);
  }
  if (delta <= 1.5e0 * rho) delta = rho;

  // Set KNEW to the index of the next interpolation point to be deleted.
  rhosq = std::max(tenth * delta, rho) * std::max(tenth * delta, rho);
  ktemp = 0;
  detrat = zero;
  if (f >= fsave) {
    ktemp = kopt;
    detrat = one;
  }
  for (k = 1; k <= npt; k++) {
    hdiag = zero;
    for (j = 1; j <= nptm; j++) {
      temp = one;
      if (j < idz) temp = -one;
      hdiag = hdiag + temp * ZMAT(k, j) * ZMAT(k, j);
    }
    temp = std::fabs(beta * hdiag + VLAG(k) * VLAG(k));
    distsq = zero;
    for (j = 1; j <= n; j++) {
      distsq = distsq + (XPT(k, j) - XOPT(j)) * (XPT(k, j) - XOPT(j));
    }
    if (distsq > rhosq) {
      temp = temp * ((distsq / rhosq) * (distsq / rhosq) * (distsq / rhosq));
    }
    if (temp > detrat && k != ktemp) {
      detrat = temp;
      knew = k;
    }
  }
  if (knew == 0) goto L460;

  // Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
  // can be moved. Begin the updating of the quadratic model, starting
  // with the explicit second derivative term.
L410:
  update(n, npt, bmat, zmat, idz, ndim, vlag, beta, knew, w);
  FVAL(knew) = f;
  ih = 0;
  for (i = 1; i <= n; i++) {
    temp = PQ(knew) * XPT(knew, i);
    for (j = 1; j <= i; j++) {
      ih = ih + 1;
      HQ(ih) = HQ(ih) + temp * XPT(knew, j);
    }
  }
  PQ(knew) = zero;

  // Update the other second derivative parameters, and then the gradient
  // vector of the model. Also include the new interpolation point.
  for (j = 1; j <= nptm; j++) {
    temp = diff * ZMAT(knew, j);
    if (j < idz) temp = -temp;
    for (k = 1; k <= npt; k++) {
      PQ(k) = PQ(k) + temp * ZMAT(k, j);
    }
  }
  gqsq = zero;
  for (i = 1; i <= n; i++) {
    GQ(i) = GQ(i) + diff * BMAT(knew, i);
    gqsq = gqsq + GQ(i) * GQ(i);
    XPT(knew, i) = XNEW(i);
  }

  // If a trust region step makes a small change to the objective function,
  // then calculate the gradient of the least Frobenius norm interpolant at
  // XBASE, and store it in W, using VLAG for a vector of right hand sides.
  if (ksave == 0 && delta == rho) {
    if (std::fabs(ratio) > 1.0e-2) {
      itest = 0;
    } else {
      for (k = 1; k <= npt; k++) {
        VLAG(k) = FVAL(k) - FVAL(kopt);
      }
      gisq = zero;
      for (i = 1; i <= n; i++) {
        sum = zero;
        for (k = 1; k <= npt; k++) {
          sum = sum + BMAT(k, i) * VLAG(k);
        }
        gisq = gisq + sum * sum;
        W(i) = sum;
      }

      // Test whether to replace the new quadratic model by the least
      // Frobenius norm interpolant, making the replacement if the test is
      // satisfied.
      itest = itest + 1;
      if (gqsq < 1.0e2 * gisq) itest = 0;
      if (itest >= 3) {
        for (i = 1; i <= n; i++) {
          GQ(i) = W(i);
        }
        for (ih = 1; ih <= nh; ih++) {
          HQ(ih) = zero;
        }
        for (j = 1; j <= nptm; j++) {
          W(j) = zero;
          for (k = 1; k <= npt; k++) {
            W(j) = W(j) + VLAG(k) * ZMAT(k, j);
          }
          if (j < idz) W(j) = -W(j);
        }
        for (k = 1; k <= npt; k++) {
          PQ(k) = zero;
          for (j = 1; j <= nptm; j++) {
            PQ(k) = PQ(k) + ZMAT(k, j) * W(j);
          }
        }
        itest = 0;
      }
    }
  }
  if (f < fsave) kopt = knew;

  // If a trust region step has provided a sufficient decrease in F, then
  // branch for another trust region calculation. The case KSAVE>0 occurs
  // when the new function value was calculated by a model step.
  if (f <= fsave + tenth * vquad) goto L100;
  if (ksave > 0) goto L100;

  // Alternatively, find out if the interpolation points are close enough
  // to the best point so far.
  knew = 0;
L460:
  distsq = 4.0e0 * delta * delta;
  for (k = 1; k <= npt; k++) {
    sum = zero;
    for (j = 1; j <= n; j++) {
      sum = sum + (XPT(k, j) - XOPT(j)) * (XPT(k, j) - XOPT(j));
    }
    if (sum > distsq) {
      knew = k;
      distsq = sum;
    }
  }

  // If KNEW is positive, then set DSTEP, and branch back for the next
  // iteration, which will generate a "model step".
  if (knew > 0) {
    dstep = std::max(std::min(tenth * std::sqrt(distsq), half * delta), rho);
    dsq = dstep * dstep;
    goto L120;
  }
  if (ratio > zero) goto L100;
  if (std::max(delta, dnorm) > rho) goto L100;

  // The calculations with the current value of RHO are complete. Pick the
  // next values of RHO and DELTA.
L490:
  if (rho > rhoend) {
    delta = half * rho;
    ratio = rho / rhoend;
    if (ratio <= 16.0e0) {
      rho = rhoend;
    } else if (ratio <= 250.0e0) {
      rho = std::sqrt(ratio) * rhoend;
    } else {
      rho = tenth * rho;
    }
    delta = std::max(delta, rho);
    if (iprint >= 2) {
      calfun.minqit(rho, nf, fopt, xbase, xopt);
    }
    goto L90;
  }

  // Return from the calculation, after another Newton-Raphson step, if
  // it is too short to have been tried before.
  if (knew == -1) goto L290;
L530:
  if (fopt <= f) {
    for (i = 1; i <= n; i++) {
      X(i) = XBASE(i) + XOPT(i);
    }
    f = fopt;
  }
  if (iprint >= 1) {
    calfun.minqir(f, nf, x);
  }
  return;
#undef X
#undef XBASE
#undef XOPT
#undef XNEW
#undef XPT
#undef FVAL
#undef GQ
#undef HQ
#undef PQ
#undef BMAT
#undef ZMAT
#undef D
#undef VLAG
#undef W
}

// ---------------------------------------------------------------------------
// newuoa.f: SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,IERR)
// ---------------------------------------------------------------------------
void newuoa(int n, int npt, double *x, double rhobeg, double rhoend,
            int iprint, int maxfun, double *w, int &ierr,
            MinqaCalfun &calfun) {
#define W(I) w[(I) - 1]
  int np, nptm, ndim, ixb, ixo, ixn, ixp, ifv, igq, ihq, ipq, ibmat, izmat,
      id, ivl, iw;

  // Partition the working space array, so that different parts of it can be
  // treated separately by the subroutine that performs the main calculation.
  np = n + 1;
  nptm = npt - np;
  if (npt < n + 2 || npt > ((n + 2) * np) / 2) {
    ierr = 10;
    goto L20;
  }
  ndim = npt + n;
  ixb = 1;
  ixo = ixb + n;
  ixn = ixo + n;
  ixp = ixn + n;
  ifv = ixp + n * npt;
  igq = ifv + npt;
  ihq = igq + n;
  ipq = ihq + (n * np) / 2;
  ibmat = ipq + npt;
  izmat = ibmat + ndim * n;
  id = izmat + npt * nptm;
  ivl = id + n;
  iw = ivl + ndim;

  // The above settings provide a partition of W for subroutine NEWUOB.
  // The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
  // W plus the space that is needed by the last array of NEWUOB.
  newuob(n, npt, x, rhobeg, rhoend, iprint, maxfun, &W(ixb), &W(ixo),
         &W(ixn), &W(ixp), &W(ifv), &W(igq), &W(ihq), &W(ipq), &W(ibmat),
         &W(izmat), ndim, &W(id), &W(ivl), &W(iw), ierr, calfun);
L20:
  return;
#undef W
}

}  // namespace

// Thread-safe replacement for minqa.cpp's newuoa_cpp().
extern "C" int newuoa_solve_c(int n, const double *par,
                              minqa_c_objfun_t objfun, void *userdata,
                              const minqa_options_t *opts,
                              minqa_result_t *result) {
  minqa_options_t defopts;
  if (opts == nullptr && n > 0 && par != nullptr) {
    defopts = minqa_options_default(n, par);
    opts = &defopts;
  }
  if (opts == nullptr) {
    if (result != nullptr) {
      minqa_result_zero(result);
      result->n = n;
      result->ierr = MINQA_ERR_INVALID;
    }
    return MINQA_ERR_INVALID;
  }
  int chk = minqa_check_args(n, par, objfun, opts, result);
  if (chk != MINQA_OK) return chk;
  MinqaCalfun calfun(objfun, userdata, n, opts);
  std::vector<double> x, w;
  int np = opts->npt;
  try {
    x.assign(par, par + n);
    // Same size as minqa.cpp: (np+13)*(np+n)+(3*n*(n+3))/2. It is only
    // clamped to at least one element so that an invalid (negative) npt
    // still reaches NEWUOA's IERR = 10 check instead of a length_error.
    long long wlen = (static_cast<long long>(np) + 13) *
                         (static_cast<long long>(np) + n) +
                     (3LL * n * (n + 3)) / 2;
    w.assign(static_cast<size_t>(std::max(wlen, 1LL)), 0.0);
  } catch (const std::bad_alloc &) {
    minqa_result_zero(result);
    result->n = n;
    result->ierr = MINQA_ERR_NOMEM;
    return MINQA_ERR_NOMEM;
  }
  return minqa_run(n, x.data(), calfun, result, [&]() -> int {
    int ierr = 0;
    newuoa(n, np, x.data(), opts->rhobeg, opts->rhoend, opts->iprint,
           opts->maxfun, w.data(), ierr, calfun);
    return ierr;
  });
}
