// Thread-safe C++17 port of M. J. D. Powell's BOBYQA as shipped (with
// R-specific edits by John Nash et al.) in the CRAN package minqa 1.2.8:
// bobyqa.f, bobyqb.f, prelim.f, rescue.f, trsbox.f, altmov.f and
// updatebobyqa.f. The translation is line by line: Powell's variable names,
// the partitioning of the work array W, column-major 2-D arrays accessed
// with 1-based index macros, labels as gotos and the Fortran order of
// floating point operations are all kept, so results are bitwise identical
// to the Fortran compiled without FMA/fast-math.
//
// Differences from the Fortran:
//  * CALFUN(N,X,IPRINT) is MinqaCalfun::operator()(x); minqit/minqir are
//    MinqaCalfun::minqit/minqir (printing via an optional callback).
//  * RESCUE in minqa calls CALFUN(N,X,IPRINT) where X is an undeclared
//    (implicitly scalar, uninitialized) local, i.e. undefined behaviour.
//    Powell's original code evaluates F at W(1..N), which RESCUE has just
//    filled with the new point; this port does the same.
//  * Locals that the Fortran may read before assignment on paths that are
//    unreachable in practice are zero-initialized (KSAV in ALTMOV is
//    initialized to KOPT so that XNEW=XOPT in that degenerate case).
//
// No static or global mutable state, no R API: safe to call concurrently.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <new>
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

// Hook used only by the validation harness to count RESCUE calls.
#ifndef MINQA_BOBYQA_RESCUE_HOOK
#define MINQA_BOBYQA_RESCUE_HOOK()
#endif

// 1-based, column-major array accessors mirroring the Fortran declarations.
#define X(i) x[(i) - 1]
#define XL(i) xl[(i) - 1]
#define XU(i) xu[(i) - 1]
#define W(i) w[(i) - 1]
#define XBASE(i) xbase[(i) - 1]
#define XPT(k, j) xpt[((k) - 1) + ((j) - 1) * npt]
#define FVAL(i) fval[(i) - 1]
#define XOPT(i) xopt[(i) - 1]
#define GOPT(i) gopt[(i) - 1]
#define HQ(i) hq[(i) - 1]
#define PQ(i) pq[(i) - 1]
#define BMAT(i, j) bmat[((i) - 1) + ((j) - 1) * ndim]
#define ZMAT(k, j) zmat[((k) - 1) + ((j) - 1) * npt]
#define SL(i) sl[(i) - 1]
#define SU(i) su[(i) - 1]
#define XNEW(i) xnew[(i) - 1]
#define XALT(i) xalt[(i) - 1]
#define D(i) d[(i) - 1]
#define VLAG(i) vlag[(i) - 1]
#define GNEW(i) gnew[(i) - 1]
#define XBDI(i) xbdi[(i) - 1]
#define S(i) s[(i) - 1]
#define HS(i) hs[(i) - 1]
#define HRED(i) hred[(i) - 1]
#define GLAG(i) glag[(i) - 1]
#define HCOL(i) hcol[(i) - 1]
#define PTSAUX(i, j) ptsaux[((i) - 1) + ((j) - 1) * 2]
#define PTSID(i) ptsid[(i) - 1]

namespace {

// ---------------------------------------------------------------------------
// updatebobyqa.f: SUBROUTINE UPDATEBOBYQA
// ---------------------------------------------------------------------------
void updatebobyqa(int n, int npt, double *bmat, double *zmat, int ndim,
                  double *vlag, double beta, double denom, int knew,
                  double *w) {
  double one, zero, ztest, temp, tempa, tempb, alpha, tau;
  int nptm, k, j, i, jp;

  one = 1.0;
  zero = 0.0;
  nptm = npt - n - 1;
  ztest = zero;
  for (k = 1; k <= npt; k++) {
    for (j = 1; j <= nptm; j++) {
      ztest = std::max(ztest, std::fabs(ZMAT(k, j)));
    }
  }
  ztest = 1.0e-20 * ztest;
  //
  //     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
  //
  for (j = 2; j <= nptm; j++) {
    if (std::fabs(ZMAT(knew, j)) > ztest) {
      temp = std::sqrt(ZMAT(knew, 1) * ZMAT(knew, 1) +
                       ZMAT(knew, j) * ZMAT(knew, j));
      tempa = ZMAT(knew, 1) / temp;
      tempb = ZMAT(knew, j) / temp;
      for (i = 1; i <= npt; i++) {
        temp = tempa * ZMAT(i, 1) + tempb * ZMAT(i, j);
        ZMAT(i, j) = tempa * ZMAT(i, j) - tempb * ZMAT(i, 1);
        ZMAT(i, 1) = temp;
      }
    }
    ZMAT(knew, j) = zero;
  }
  //
  //     Put the first NPT components of the KNEW-th column of HLAG into W,
  //     and calculate the parameters of the updating formula.
  //
  for (i = 1; i <= npt; i++) {
    W(i) = ZMAT(knew, 1) * ZMAT(i, 1);
  }
  alpha = W(knew);
  tau = VLAG(knew);
  VLAG(knew) = VLAG(knew) - one;
  //
  //     Complete the updating of ZMAT.
  //
  temp = std::sqrt(denom);
  tempb = ZMAT(knew, 1) / temp;
  tempa = tau / temp;
  for (i = 1; i <= npt; i++) {
    ZMAT(i, 1) = tempa * ZMAT(i, 1) - tempb * VLAG(i);
  }
  //
  //     Finally, update the matrix BMAT.
  //
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
}

// ---------------------------------------------------------------------------
// prelim.f: SUBROUTINE PRELIM
// ---------------------------------------------------------------------------
void prelim(int n, int npt, double *x, const double *xl, const double *xu,
            double rhobeg, int maxfun, double *xbase, double *xpt,
            double *fval, double *gopt, double *hq, double *pq, double *bmat,
            double *zmat, int ndim, double *sl, double *su, int &nf,
            int &kopt, MinqaCalfun &calfun) {
  double half, one, two, zero, rhosq, recip, stepa = 0.0, stepb = 0.0, f,
                                                 fbeg = 0.0, temp, diff;
  int np, j, k, i, ih, nfm, nfx, itemp, jpt = 0, ipt = 0;

  half = 0.5;
  one = 1.0;
  two = 2.0;
  zero = 0.0;
  rhosq = rhobeg * rhobeg;
  recip = one / rhosq;
  np = n + 1;
  //
  //     Set XBASE to the initial vector of variables, and set the initial
  //     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
  //
  for (j = 1; j <= n; j++) {
    XBASE(j) = X(j);
    for (k = 1; k <= npt; k++) {
      XPT(k, j) = zero;
    }
    for (i = 1; i <= ndim; i++) {
      BMAT(i, j) = zero;
    }
  }
  for (ih = 1; ih <= (n * np) / 2; ih++) {
    HQ(ih) = zero;
  }
  for (k = 1; k <= npt; k++) {
    PQ(k) = zero;
    for (j = 1; j <= npt - np; j++) {
      ZMAT(k, j) = zero;
    }
  }
  //
  //     Begin the initialization procedure. NF becomes one more than the
  //     number of function values so far. The coordinates of the displacement
  //     of the next initial interpolation point from XBASE are set in
  //     XPT(NF+1,.).
  //
  nf = 0;
L50:
  nfm = nf;
  nfx = nf - n;
  nf = nf + 1;
  if (nfm <= 2 * n) {
    if (nfm >= 1 && nfm <= n) {
      stepa = rhobeg;
      if (SU(nfm) == zero) stepa = -stepa;
      XPT(nf, nfm) = stepa;
    } else if (nfm > n) {
      stepa = XPT(nf - n, nfx);
      stepb = -rhobeg;
      if (SL(nfx) == zero) stepb = std::min(two * rhobeg, SU(nfx));
      if (SU(nfx) == zero) stepb = std::max(-two * rhobeg, SL(nfx));
      XPT(nf, nfx) = stepb;
    }
  } else {
    itemp = (nfm - np) / n;
    jpt = nfm - itemp * n - n;
    ipt = jpt + itemp;
    if (ipt > n) {
      itemp = jpt;
      jpt = ipt - n;
      ipt = itemp;
    }
    XPT(nf, ipt) = XPT(ipt + 1, ipt);
    XPT(nf, jpt) = XPT(jpt + 1, jpt);
  }
  //
  //     Calculate the next value of F. The least function value so far and
  //     its index are required.
  //
  for (j = 1; j <= n; j++) {
    X(j) = std::min(std::max(XL(j), XBASE(j) + XPT(nf, j)), XU(j));
    if (XPT(nf, j) == SL(j)) X(j) = XL(j);
    if (XPT(nf, j) == SU(j)) X(j) = XU(j);
  }
  f = calfun(x);
  FVAL(nf) = f;
  if (nf == 1) {
    fbeg = f;
    kopt = 1;
  } else if (f < FVAL(kopt)) {
    kopt = nf;
  }
  //
  //     Set the nonzero initial elements of BMAT and the quadratic model in
  //     the cases when NF is at most 2*N+1. If NF exceeds N+1, then the
  //     positions of the NF-th and (NF-N)-th interpolation points may be
  //     switched, in order that the function value at the first of them
  //     contributes to the off-diagonal second derivative terms of the
  //     initial quadratic model.
  //
  if (nf <= 2 * n + 1) {
    if (nf >= 2 && nf <= n + 1) {
      GOPT(nfm) = (f - fbeg) / stepa;
      if (npt < nf + n) {
        BMAT(1, nfm) = -one / stepa;
        BMAT(nf, nfm) = one / stepa;
        BMAT(npt + nfm, nfm) = -half * rhosq;
      }
    } else if (nf >= n + 2) {
      ih = (nfx * (nfx + 1)) / 2;
      temp = (f - fbeg) / stepb;
      diff = stepb - stepa;
      HQ(ih) = two * (temp - GOPT(nfx)) / diff;
      GOPT(nfx) = (GOPT(nfx) * stepb - temp * stepa) / diff;
      if (stepa * stepb < zero) {
        if (f < FVAL(nf - n)) {
          FVAL(nf) = FVAL(nf - n);
          FVAL(nf - n) = f;
          if (kopt == nf) kopt = nf - n;
          XPT(nf - n, nfx) = stepb;
          XPT(nf, nfx) = stepa;
        }
      }
      BMAT(1, nfx) = -(stepa + stepb) / (stepa * stepb);
      BMAT(nf, nfx) = -half / XPT(nf - n, nfx);
      BMAT(nf - n, nfx) = -BMAT(1, nfx) - BMAT(nf, nfx);
      ZMAT(1, nfx) = std::sqrt(two) / (stepa * stepb);
      ZMAT(nf, nfx) = std::sqrt(half) / rhosq;
      ZMAT(nf - n, nfx) = -ZMAT(1, nfx) - ZMAT(nf, nfx);
    }
    //
    //     Set the off-diagonal second derivatives of the Lagrange functions
    //     and the initial quadratic model.
    //
  } else {
    ih = (ipt * (ipt - 1)) / 2 + jpt;
    ZMAT(1, nfx) = recip;
    ZMAT(nf, nfx) = recip;
    ZMAT(ipt + 1, nfx) = -recip;
    ZMAT(jpt + 1, nfx) = -recip;
    temp = XPT(nf, ipt) * XPT(nf, jpt);
    HQ(ih) = (fbeg - FVAL(ipt + 1) - FVAL(jpt + 1) + f) / temp;
  }
  if (nf < npt && nf < maxfun) goto L50;
}

// ---------------------------------------------------------------------------
// trsbox.f: SUBROUTINE TRSBOX
// ---------------------------------------------------------------------------
void trsbox(int n, int npt, double *xpt, double *xopt, double *gopt,
            double *hq, double *pq, double *sl, double *su, double delta,
            double *xnew, double *d, double *gnew, double *xbdi, double *s,
            double *hs, double *hred, double &dsq, double &crvmin) {
  double half, one, onemin, zero, delsq, qred, beta, stepsq,
      gredsq = 0.0, resid, ds, shs, temp, blen, stplen, xsum, sdec,
      ggsav = 0.0, dredsq = 0.0, dredg = 0.0, sredg = 0.0, angbd = 0.0,
      tempa, tempb, ssq, xsav = 0.0, dhs, dhd, redmax, redsav, angt = 0.0,
      sth, rednew, rdprev = 0.0, rdnext = 0.0, cth;
  int iterc, nact, i, itermax = 0, iact = 0, itcsav = 0, isav, iu, ih, j, k;

  half = 0.5;
  one = 1.0;
  onemin = -1.0;
  zero = 0.0;
  //
  //     The sign of GOPT(I) gives the sign of the change to the I-th variable
  //     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
  //     or not to fix the I-th variable at one of its bounds initially, with
  //     NACT being set to the number of fixed variables. D and GNEW are also
  //     set for the first iteration. DELSQ is the upper bound on the sum of
  //     squares of the free variables. QRED is the reduction in Q so far.
  //
  iterc = 0;
  nact = 0;
  for (i = 1; i <= n; i++) {
    XBDI(i) = zero;
    if (XOPT(i) <= SL(i)) {
      if (GOPT(i) >= zero) XBDI(i) = onemin;
    } else if (XOPT(i) >= SU(i)) {
      if (GOPT(i) <= zero) XBDI(i) = one;
    }
    if (XBDI(i) != zero) nact = nact + 1;
    D(i) = zero;
    GNEW(i) = GOPT(i);
  }
  delsq = delta * delta;
  qred = zero;
  crvmin = onemin;
  //
  //     Set the next search direction of the conjugate gradient method. It is
  //     the steepest descent direction initially and when the iterations are
  //     restarted because a variable has just been fixed by a bound, and of
  //     course the components of the fixed variables are zero. ITERMAX is an
  //     upper bound on the indices of the conjugate gradient iterations.
  //
L20:
  beta = zero;
L30:
  stepsq = zero;
  for (i = 1; i <= n; i++) {
    if (XBDI(i) != zero) {
      S(i) = zero;
    } else if (beta == zero) {
      S(i) = -GNEW(i);
    } else {
      S(i) = beta * S(i) - GNEW(i);
    }
    stepsq = stepsq + S(i) * S(i);
  }
  if (stepsq == zero) goto L190;
  if (beta == zero) {
    gredsq = stepsq;
    itermax = iterc + n - nact;
  }
  if (gredsq * delsq <= 1.0e-4 * qred * qred) goto L190;
  //
  //     Multiply the search direction by the second derivative matrix of Q
  //     and calculate some scalars for the choice of steplength. Then set
  //     BLEN to the length of the the step to the trust region boundary and
  //     STPLEN to the steplength, ignoring the simple bounds.
  //
  goto L210;
L50:
  resid = delsq;
  ds = zero;
  shs = zero;
  for (i = 1; i <= n; i++) {
    if (XBDI(i) == zero) {
      resid = resid - D(i) * D(i);
      ds = ds + S(i) * D(i);
      shs = shs + S(i) * HS(i);
    }
  }
  if (resid <= zero) goto L90;
  temp = std::sqrt(stepsq * resid + ds * ds);
  if (ds < zero) {
    blen = (temp - ds) / stepsq;
  } else {
    blen = resid / (temp + ds);
  }
  stplen = blen;
  if (shs > zero) {
    stplen = std::min(blen, gredsq / shs);
  }
  //
  //     Reduce STPLEN if necessary in order to preserve the simple bounds,
  //     letting IACT be the index of the new constrained variable.
  //
  iact = 0;
  for (i = 1; i <= n; i++) {
    if (S(i) != zero) {
      xsum = XOPT(i) + D(i);
      if (S(i) > zero) {
        temp = (SU(i) - xsum) / S(i);
      } else {
        temp = (SL(i) - xsum) / S(i);
      }
      if (temp < stplen) {
        stplen = temp;
        iact = i;
      }
    }
  }
  //
  //     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
  //
  sdec = zero;
  if (stplen > zero) {
    iterc = iterc + 1;
    temp = shs / stepsq;
    if (iact == 0 && temp > zero) {
      crvmin = std::min(crvmin, temp);
      if (crvmin == onemin) crvmin = temp;
    }
    ggsav = gredsq;
    gredsq = zero;
    for (i = 1; i <= n; i++) {
      GNEW(i) = GNEW(i) + stplen * HS(i);
      if (XBDI(i) == zero) gredsq = gredsq + GNEW(i) * GNEW(i);
      D(i) = D(i) + stplen * S(i);
    }
    sdec = std::max(stplen * (ggsav - half * stplen * shs), zero);
    qred = qred + sdec;
  }
  //
  //     Restart the conjugate gradient method if it has hit a new bound.
  //
  if (iact > 0) {
    nact = nact + 1;
    XBDI(iact) = one;
    if (S(iact) < zero) XBDI(iact) = onemin;
    delsq = delsq - D(iact) * D(iact);
    if (delsq <= zero) goto L90;
    goto L20;
  }
  //
  //     If STPLEN is less than BLEN, then either apply another conjugate
  //     gradient iteration or RETURN.
  //
  if (stplen < blen) {
    if (iterc == itermax) goto L190;
    if (sdec <= 0.01 * qred) goto L190;
    beta = gredsq / ggsav;
    goto L30;
  }
L90:
  crvmin = zero;
  //
  //     Prepare for the alternative iteration by calculating some scalars
  //     and by multiplying the reduced D by the second derivative matrix of
  //     Q, where S holds the reduced D in the call of GGMULT.
  //
L100:
  if (nact >= n - 1) goto L190;
  dredsq = zero;
  dredg = zero;
  gredsq = zero;
  for (i = 1; i <= n; i++) {
    if (XBDI(i) == zero) {
      dredsq = dredsq + D(i) * D(i);
      dredg = dredg + D(i) * GNEW(i);
      gredsq = gredsq + GNEW(i) * GNEW(i);
      S(i) = D(i);
    } else {
      S(i) = zero;
    }
  }
  itcsav = iterc;
  goto L210;
  //
  //     Let the search direction S be a linear combination of the reduced D
  //     and the reduced G that is orthogonal to the reduced D.
  //
L120:
  iterc = iterc + 1;
  temp = gredsq * dredsq - dredg * dredg;
  if (temp <= 1.0e-4 * qred * qred) goto L190;
  temp = std::sqrt(temp);
  for (i = 1; i <= n; i++) {
    if (XBDI(i) == zero) {
      S(i) = (dredg * D(i) - dredsq * GNEW(i)) / temp;
    } else {
      S(i) = zero;
    }
  }
  sredg = -temp;
  //
  //     By considering the simple bounds on the variables, calculate an upper
  //     bound on the tangent of half the angle of the alternative iteration,
  //     namely ANGBD, except that, if already a free variable has reached a
  //     bound, there is a branch back to label 100 after fixing that
  //     variable.
  //
  angbd = one;
  iact = 0;
  for (i = 1; i <= n; i++) {
    if (XBDI(i) == zero) {
      tempa = XOPT(i) + D(i) - SL(i);
      tempb = SU(i) - XOPT(i) - D(i);
      if (tempa <= zero) {
        nact = nact + 1;
        XBDI(i) = onemin;
        goto L100;
      } else if (tempb <= zero) {
        nact = nact + 1;
        XBDI(i) = one;
        goto L100;
      }
      ssq = D(i) * D(i) + S(i) * S(i);
      temp = ssq - (XOPT(i) - SL(i)) * (XOPT(i) - SL(i));
      if (temp > zero) {
        temp = std::sqrt(temp) - S(i);
        if (angbd * temp > tempa) {
          angbd = tempa / temp;
          iact = i;
          xsav = onemin;
        }
      }
      temp = ssq - (SU(i) - XOPT(i)) * (SU(i) - XOPT(i));
      if (temp > zero) {
        temp = std::sqrt(temp) + S(i);
        if (angbd * temp > tempb) {
          angbd = tempb / temp;
          iact = i;
          xsav = one;
        }
      }
    }
  }
  //
  //     Calculate HHD and some curvatures for the alternative iteration.
  //
  goto L210;
L150:
  shs = zero;
  dhs = zero;
  dhd = zero;
  for (i = 1; i <= n; i++) {
    if (XBDI(i) == zero) {
      shs = shs + S(i) * HS(i);
      dhs = dhs + D(i) * HS(i);
      dhd = dhd + D(i) * HRED(i);
    }
  }
  //
  //     Seek the greatest reduction in Q for a range of equally spaced values
  //     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
  //     the alternative iteration.
  //
  redmax = zero;
  isav = 0;
  redsav = zero;
  iu = static_cast<int>(17.0 * angbd + 3.1);
  for (i = 1; i <= iu; i++) {
    angt = angbd * static_cast<double>(i) / static_cast<double>(iu);
    sth = (angt + angt) / (one + angt * angt);
    temp = shs + angt * (angt * dhd - dhs - dhs);
    rednew = sth * (angt * dredg - sredg - half * sth * temp);
    if (rednew > redmax) {
      redmax = rednew;
      isav = i;
      rdprev = redsav;
    } else if (i == isav + 1) {
      rdnext = rednew;
    }
    redsav = rednew;
  }
  //
  //     Return if the reduction is zero. Otherwise, set the sine and cosine
  //     of the angle of the alternative iteration, and calculate SDEC.
  //
  if (isav == 0) goto L190;
  if (isav < iu) {
    temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
    angt = angbd * (static_cast<double>(isav) + half * temp) /
           static_cast<double>(iu);
  }
  cth = (one - angt * angt) / (one + angt * angt);
  sth = (angt + angt) / (one + angt * angt);
  temp = shs + angt * (angt * dhd - dhs - dhs);
  sdec = sth * (angt * dredg - sredg - half * sth * temp);
  if (sdec <= zero) goto L190;
  //
  //     Update GNEW, D and HRED. If the angle of the alternative iteration
  //     is restricted by a bound on a free variable, that variable is fixed
  //     at the bound.
  //
  dredg = zero;
  gredsq = zero;
  for (i = 1; i <= n; i++) {
    GNEW(i) = GNEW(i) + (cth - one) * HRED(i) + sth * HS(i);
    if (XBDI(i) == zero) {
      D(i) = cth * D(i) + sth * S(i);
      dredg = dredg + D(i) * GNEW(i);
      gredsq = gredsq + GNEW(i) * GNEW(i);
    }
    HRED(i) = cth * HRED(i) + sth * HS(i);
  }
  qred = qred + sdec;
  if (iact > 0 && isav == iu) {
    nact = nact + 1;
    XBDI(iact) = xsav;
    goto L100;
  }
  //
  //     If SDEC is sufficiently small, then RETURN after setting XNEW to
  //     XOPT+D, giving careful attention to the bounds.
  //
  if (sdec > 0.01 * qred) goto L120;
L190:
  dsq = zero;
  for (i = 1; i <= n; i++) {
    XNEW(i) = std::max(std::min(XOPT(i) + D(i), SU(i)), SL(i));
    if (XBDI(i) == onemin) XNEW(i) = SL(i);
    if (XBDI(i) == one) XNEW(i) = SU(i);
    D(i) = XNEW(i) - XOPT(i);
    dsq = dsq + D(i) * D(i);
  }
  return;

  //     The following instructions multiply the current S-vector by the
  //     second derivative matrix of the quadratic model, putting the product
  //     in HS. They are reached from three different parts of the software
  //     above and they can be regarded as an external subroutine.
  //
L210:
  ih = 0;
  for (j = 1; j <= n; j++) {
    HS(j) = zero;
    for (i = 1; i <= j; i++) {
      ih = ih + 1;
      if (i < j) HS(j) = HS(j) + HQ(ih) * S(i);
      HS(i) = HS(i) + HQ(ih) * S(j);
    }
  }
  for (k = 1; k <= npt; k++) {
    if (PQ(k) != zero) {
      temp = zero;
      for (j = 1; j <= n; j++) {
        temp = temp + XPT(k, j) * S(j);
      }
      temp = temp * PQ(k);
      for (i = 1; i <= n; i++) {
        HS(i) = HS(i) + temp * XPT(k, i);
      }
    }
  }
  if (crvmin != zero) goto L50;
  if (iterc > itcsav) goto L150;
  for (i = 1; i <= n; i++) {
    HRED(i) = HS(i);
  }
  goto L120;
}

// ---------------------------------------------------------------------------
// altmov.f: SUBROUTINE ALTMOV
// ---------------------------------------------------------------------------
void altmov(int n, int npt, double *xpt, double *xopt, double *bmat,
            double *zmat, int ndim, double *sl, double *su, int kopt,
            int knew, double adelt, double *xnew, double *xalt,
            double &alpha, double &cauchy, double *glag, double *hcol,
            double *w) {
  double half, one, zero, cons, temp, ha, presav, dderiv, distsq, subd,
      slbd, sumin, diff, step = 0.0, vlag, tempd, tempa, tempb, predsq,
      stpsav = 0.0, bigstp, wfixsq, ggfree, wsqsav, gw, curv, scale,
      csave = 0.0;
  int k, j, i, ilbd, iubd, isbd, ksav, ibdsav, iflag;

  half = 0.5;
  one = 1.0;
  zero = 0.0;
  cons = one + std::sqrt(2.0);
  ksav = kopt;
  for (k = 1; k <= npt; k++) {
    HCOL(k) = zero;
  }
  for (j = 1; j <= npt - n - 1; j++) {
    temp = ZMAT(knew, j);
    for (k = 1; k <= npt; k++) {
      HCOL(k) = HCOL(k) + temp * ZMAT(k, j);
    }
  }
  alpha = HCOL(knew);
  ha = half * alpha;
  //
  //     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
  //
  for (i = 1; i <= n; i++) {
    GLAG(i) = BMAT(knew, i);
  }
  for (k = 1; k <= npt; k++) {
    temp = zero;
    for (j = 1; j <= n; j++) {
      temp = temp + XPT(k, j) * XOPT(j);
    }
    temp = HCOL(k) * temp;
    for (i = 1; i <= n; i++) {
      GLAG(i) = GLAG(i) + temp * XPT(k, i);
    }
  }
  //
  //     Search for a large denominator along the straight lines through XOPT
  //     and another interpolation point. SLBD and SUBD will be lower and
  //     upper bounds on the step along each of these lines in turn. PREDSQ
  //     will be set to the square of the predicted denominator for each line.
  //     PRESAV will be set to the largest admissible value of PREDSQ that
  //     occurs.
  //
  presav = zero;
  for (k = 1; k <= npt; k++) {
    if (k == kopt) continue;  // GOTO 80
    dderiv = zero;
    distsq = zero;
    for (i = 1; i <= n; i++) {
      temp = XPT(k, i) - XOPT(i);
      dderiv = dderiv + GLAG(i) * temp;
      distsq = distsq + temp * temp;
    }
    subd = adelt / std::sqrt(distsq);
    slbd = -subd;
    ilbd = 0;
    iubd = 0;
    sumin = std::min(one, subd);
    //
    //     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
    //
    for (i = 1; i <= n; i++) {
      temp = XPT(k, i) - XOPT(i);
      if (temp > zero) {
        if (slbd * temp < SL(i) - XOPT(i)) {
          slbd = (SL(i) - XOPT(i)) / temp;
          ilbd = -i;
        }
        if (subd * temp > SU(i) - XOPT(i)) {
          subd = std::max(sumin, (SU(i) - XOPT(i)) / temp);
          iubd = i;
        }
      } else if (temp < zero) {
        if (slbd * temp > SU(i) - XOPT(i)) {
          slbd = (SU(i) - XOPT(i)) / temp;
          ilbd = i;
        }
        if (subd * temp < SL(i) - XOPT(i)) {
          subd = std::max(sumin, (SL(i) - XOPT(i)) / temp);
          iubd = -i;
        }
      }
    }
    //
    //     Seek a large modulus of the KNEW-th Lagrange function when the index
    //     of the other interpolation point on the line through XOPT is KNEW.
    //
    if (k == knew) {
      diff = dderiv - one;
      step = slbd;
      vlag = slbd * (dderiv - slbd * diff);
      isbd = ilbd;
      temp = subd * (dderiv - subd * diff);
      if (std::fabs(temp) > std::fabs(vlag)) {
        step = subd;
        vlag = temp;
        isbd = iubd;
      }
      tempd = half * dderiv;
      tempa = tempd - diff * slbd;
      tempb = tempd - diff * subd;
      if (tempa * tempb < zero) {
        temp = tempd * tempd / diff;
        if (std::fabs(temp) > std::fabs(vlag)) {
          step = tempd / diff;
          vlag = temp;
          isbd = 0;
        }
      }
      //
      //     Search along each of the other lines through XOPT and another
      //     point.
      //
    } else {
      step = slbd;
      vlag = slbd * (one - slbd);
      isbd = ilbd;
      temp = subd * (one - subd);
      if (std::fabs(temp) > std::fabs(vlag)) {
        step = subd;
        vlag = temp;
        isbd = iubd;
      }
      if (subd > half) {
        if (std::fabs(vlag) < 0.25) {
          step = half;
          vlag = 0.25;
          isbd = 0;
        }
      }
      vlag = vlag * dderiv;
    }
    //
    //     Calculate PREDSQ for the current line search and maintain PRESAV.
    //
    temp = step * (one - step) * distsq;
    predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
    if (predsq > presav) {
      presav = predsq;
      ksav = k;
      stpsav = step;
      ibdsav = isbd;
    }
  }
  //
  //     Construct XNEW in a way that satisfies the bound constraints exactly.
  //
  ibdsav = 0;
  for (i = 1; i <= n; i++) {
    temp = XOPT(i) + stpsav * (XPT(ksav, i) - XOPT(i));
    XNEW(i) = std::max(SL(i), std::min(SU(i), temp));
  }
  if (ibdsav < 0) XNEW(-ibdsav) = SL(-ibdsav);
  if (ibdsav > 0) XNEW(ibdsav) = SU(ibdsav);
  //
  //     Prepare for the iterative method that assembles the constrained
  //     Cauchy step in W. The sum of squares of the fixed components of W is
  //     formed in WFIXSQ, and the free components of W are set to BIGSTP.
  //
  bigstp = adelt + adelt;
  iflag = 0;
L100:
  wfixsq = zero;
  ggfree = zero;
  for (i = 1; i <= n; i++) {
    W(i) = zero;
    tempa = std::min(XOPT(i) - SL(i), GLAG(i));
    tempb = std::max(XOPT(i) - SU(i), GLAG(i));
    if (tempa > zero || tempb < zero) {
      W(i) = bigstp;
      ggfree = ggfree + GLAG(i) * GLAG(i);
    }
  }
  if (ggfree == zero) {
    cauchy = zero;
    return;  // GOTO 200
  }
  //
  //     Investigate whether more components of W can be fixed.
  //
L120:
  temp = adelt * adelt - wfixsq;
  if (temp > zero) {
    wsqsav = wfixsq;
    step = std::sqrt(temp / ggfree);
    ggfree = zero;
    for (i = 1; i <= n; i++) {
      if (W(i) == bigstp) {
        temp = XOPT(i) - step * GLAG(i);
        if (temp <= SL(i)) {
          W(i) = SL(i) - XOPT(i);
          wfixsq = wfixsq + W(i) * W(i);
        } else if (temp >= SU(i)) {
          W(i) = SU(i) - XOPT(i);
          wfixsq = wfixsq + W(i) * W(i);
        } else {
          ggfree = ggfree + GLAG(i) * GLAG(i);
        }
      }
    }
    if (wfixsq > wsqsav && ggfree > zero) goto L120;
  }
  //
  //     Set the remaining free components of W and all components of XALT,
  //     except that W may be scaled later.
  //
  gw = zero;
  for (i = 1; i <= n; i++) {
    if (W(i) == bigstp) {
      W(i) = -step * GLAG(i);
      XALT(i) = std::max(SL(i), std::min(SU(i), XOPT(i) + W(i)));
    } else if (W(i) == zero) {
      XALT(i) = XOPT(i);
    } else if (GLAG(i) > zero) {
      XALT(i) = SL(i);
    } else {
      XALT(i) = SU(i);
    }
    gw = gw + GLAG(i) * W(i);
  }
  //
  //     Set CURV to the curvature of the KNEW-th Lagrange function along W.
  //     Scale W by a factor less than one if that can reduce the modulus of
  //     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
  //     the square of this function.
  //
  curv = zero;
  for (k = 1; k <= npt; k++) {
    temp = zero;
    for (j = 1; j <= n; j++) {
      temp = temp + XPT(k, j) * W(j);
    }
    curv = curv + HCOL(k) * temp * temp;
  }
  if (iflag == 1) curv = -curv;
  if (curv > -gw && curv < -cons * gw) {
    scale = -gw / curv;
    for (i = 1; i <= n; i++) {
      temp = XOPT(i) + scale * W(i);
      XALT(i) = std::max(SL(i), std::min(SU(i), temp));
    }
    cauchy = (half * gw * scale) * (half * gw * scale);
  } else {
    cauchy = (gw + half * curv) * (gw + half * curv);
  }
  //
  //     If IFLAG is zero, then XALT is calculated as before after reversing
  //     the sign of GLAG. Thus two XALT vectors become available. The one
  //     that is chosen is the one that gives the larger value of CAUCHY.
  //
  if (iflag == 0) {
    for (i = 1; i <= n; i++) {
      GLAG(i) = -GLAG(i);
      W(n + i) = XALT(i);
    }
    csave = cauchy;
    iflag = 1;
    goto L100;
  }
  if (csave > cauchy) {
    for (i = 1; i <= n; i++) {
      XALT(i) = W(n + i);
    }
    cauchy = csave;
  }
}

// ---------------------------------------------------------------------------
// rescue.f: SUBROUTINE RESCUE
// ---------------------------------------------------------------------------
void rescue(int n, int npt, const double *xl, const double *xu, int maxfun,
            double *xbase, double *xpt, double *fval, double *xopt,
            double *gopt, double *hq, double *pq, double *bmat, double *zmat,
            int ndim, double *sl, double *su, int &nf, double delta,
            int &kopt, double *vlag, double *ptsaux, double *ptsid, double *w,
            MinqaCalfun &calfun) {
  double half, one, zero, sfrac, sumpq, winc, distsq, temp, fbase,
      beta = 0.0, denom = 0.0, dsqmin, sum, bsum, vlmxsq, hdiag, den,
      xp = 0.0, xq = 0.0, vquad, f, diff;
  int np, nptm, k, j, ih, i, jp, jpn, iw, ip, iq, nrem, kold, knew, kpt,
      ihp = 0, ihq;

  MINQA_BOBYQA_RESCUE_HOOK();
  half = 0.5;
  one = 1.0;
  zero = 0.0;
  np = n + 1;
  sfrac = half / static_cast<double>(np);
  nptm = npt - np;
  //
  //     Shift the interpolation points so that XOPT becomes the origin, and
  //     set the elements of ZMAT to zero. The value of SUMPQ is required in
  //     the updating of HQ below. The squares of the distances from XOPT to
  //     the other interpolation points are set at the end of W. Increments
  //     of WINC may be added later to these squares to balance the
  //     consideration of the choice of point that is going to become current.
  //
  sumpq = zero;
  winc = zero;
  for (k = 1; k <= npt; k++) {
    distsq = zero;
    for (j = 1; j <= n; j++) {
      XPT(k, j) = XPT(k, j) - XOPT(j);
      distsq = distsq + XPT(k, j) * XPT(k, j);
    }
    sumpq = sumpq + PQ(k);
    W(ndim + k) = distsq;
    winc = std::max(winc, distsq);
    for (j = 1; j <= nptm; j++) {
      ZMAT(k, j) = zero;
    }
  }
  //
  //     Update HQ so that HQ and PQ define the second derivatives of the
  //     model after XBASE has been shifted to the trust region centre.
  //
  ih = 0;
  for (j = 1; j <= n; j++) {
    W(j) = half * sumpq * XOPT(j);
    for (k = 1; k <= npt; k++) {
      W(j) = W(j) + PQ(k) * XPT(k, j);
    }
    for (i = 1; i <= j; i++) {
      ih = ih + 1;
      HQ(ih) = HQ(ih) + W(i) * XOPT(j) + W(j) * XOPT(i);
    }
  }
  //
  //     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
  //     also set the elements of PTSAUX.
  //
  for (j = 1; j <= n; j++) {
    XBASE(j) = XBASE(j) + XOPT(j);
    SL(j) = SL(j) - XOPT(j);
    SU(j) = SU(j) - XOPT(j);
    XOPT(j) = zero;
    PTSAUX(1, j) = std::min(delta, SU(j));
    PTSAUX(2, j) = std::max(-delta, SL(j));
    if (PTSAUX(1, j) + PTSAUX(2, j) < zero) {
      temp = PTSAUX(1, j);
      PTSAUX(1, j) = PTSAUX(2, j);
      PTSAUX(2, j) = temp;
    }
    if (std::fabs(PTSAUX(2, j)) < half * std::fabs(PTSAUX(1, j))) {
      PTSAUX(2, j) = half * PTSAUX(1, j);
    }
    for (i = 1; i <= ndim; i++) {
      BMAT(i, j) = zero;
    }
  }
  fbase = FVAL(kopt);
  //
  //     Set the identifiers of the artificial interpolation points that are
  //     along a coordinate direction from XOPT, and set the corresponding
  //     nonzero elements of BMAT and ZMAT.
  //
  PTSID(1) = sfrac;
  for (j = 1; j <= n; j++) {
    jp = j + 1;
    jpn = jp + n;
    PTSID(jp) = static_cast<double>(j) + sfrac;
    if (jpn <= npt) {
      PTSID(jpn) = static_cast<double>(j) / static_cast<double>(np) + sfrac;
      temp = one / (PTSAUX(1, j) - PTSAUX(2, j));
      BMAT(jp, j) = -temp + one / PTSAUX(1, j);
      BMAT(jpn, j) = temp + one / PTSAUX(2, j);
      BMAT(1, j) = -BMAT(jp, j) - BMAT(jpn, j);
      ZMAT(1, j) = std::sqrt(2.0) / std::fabs(PTSAUX(1, j) * PTSAUX(2, j));
      ZMAT(jp, j) = ZMAT(1, j) * PTSAUX(2, j) * temp;
      ZMAT(jpn, j) = -ZMAT(1, j) * PTSAUX(1, j) * temp;
    } else {
      BMAT(1, j) = -one / PTSAUX(1, j);
      BMAT(jp, j) = one / PTSAUX(1, j);
      BMAT(j + npt, j) = -half * (PTSAUX(1, j) * PTSAUX(1, j));
    }
  }
  //
  //     Set any remaining identifiers with their nonzero elements of ZMAT.
  //
  if (npt >= n + np) {
    for (k = 2 * np; k <= npt; k++) {
      iw = static_cast<int>((static_cast<double>(k - np) - half) /
                            static_cast<double>(n));
      ip = k - np - iw * n;
      iq = ip + iw;
      if (iq > n) iq = iq - n;
      PTSID(k) = static_cast<double>(ip) +
                 static_cast<double>(iq) / static_cast<double>(np) + sfrac;
      temp = one / (PTSAUX(1, ip) * PTSAUX(1, iq));
      ZMAT(1, k - np) = temp;
      ZMAT(ip + 1, k - np) = -temp;
      ZMAT(iq + 1, k - np) = -temp;
      ZMAT(k, k - np) = temp;
    }
  }
  nrem = npt;
  kold = 1;
  knew = kopt;
  //
  //     Reorder the provisional points in the way that exchanges PTSID(KOLD)
  //     with PTSID(KNEW).
  //
L80:
  for (j = 1; j <= n; j++) {
    temp = BMAT(kold, j);
    BMAT(kold, j) = BMAT(knew, j);
    BMAT(knew, j) = temp;
  }
  for (j = 1; j <= nptm; j++) {
    temp = ZMAT(kold, j);
    ZMAT(kold, j) = ZMAT(knew, j);
    ZMAT(knew, j) = temp;
  }
  PTSID(kold) = PTSID(knew);
  PTSID(knew) = zero;
  W(ndim + knew) = zero;
  nrem = nrem - 1;
  if (knew != kopt) {
    temp = VLAG(kold);
    VLAG(kold) = VLAG(knew);
    VLAG(knew) = temp;
    //
    //     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
    //     interpolation point can be changed from provisional to original.
    //     The branch to label 350 occurs if all the original points are
    //     reinstated. The nonnegative values of W(NDIM+K) are required in
    //     the search below.
    //
    updatebobyqa(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w);
    if (nrem == 0) return;  // GOTO 350
    for (k = 1; k <= npt; k++) {
      W(ndim + k) = std::fabs(W(ndim + k));
    }
  }
  //
  //     Pick the index KNEW of an original interpolation point that has not
  //     yet replaced one of the provisional interpolation points, giving
  //     attention to the closeness to XOPT and to previous tries with KNEW.
  //
L120:
  dsqmin = zero;
  for (k = 1; k <= npt; k++) {
    if (W(ndim + k) > zero) {
      if (dsqmin == zero || W(ndim + k) < dsqmin) {
        knew = k;
        dsqmin = W(ndim + k);
      }
    }
  }
  if (dsqmin == zero) goto L260;
  //
  //     Form the W-vector of the chosen original interpolation point.
  //
  for (j = 1; j <= n; j++) {
    W(npt + j) = XPT(knew, j);
  }
  for (k = 1; k <= npt; k++) {
    sum = zero;
    if (k == kopt) {
      // CONTINUE
    } else if (PTSID(k) == zero) {
      for (j = 1; j <= n; j++) {
        sum = sum + W(npt + j) * XPT(k, j);
      }
    } else {
      ip = static_cast<int>(PTSID(k));
      if (ip > 0) sum = W(npt + ip) * PTSAUX(1, ip);
      iq = static_cast<int>(static_cast<double>(np) * PTSID(k) -
                            static_cast<double>(ip * np));
      if (iq > 0) {
        iw = 1;
        if (ip == 0) iw = 2;
        sum = sum + W(npt + iq) * PTSAUX(iw, iq);
      }
    }
    W(k) = half * sum * sum;
  }
  //
  //     Calculate VLAG and BETA for the required updating of the H matrix if
  //     XPT(KNEW,.) is reinstated in the set of interpolation points.
  //
  for (k = 1; k <= npt; k++) {
    sum = zero;
    for (j = 1; j <= n; j++) {
      sum = sum + BMAT(k, j) * W(npt + j);
    }
    VLAG(k) = sum;
  }
  beta = zero;
  for (j = 1; j <= nptm; j++) {
    sum = zero;
    for (k = 1; k <= npt; k++) {
      sum = sum + ZMAT(k, j) * W(k);
    }
    beta = beta - sum * sum;
    for (k = 1; k <= npt; k++) {
      VLAG(k) = VLAG(k) + sum * ZMAT(k, j);
    }
  }
  bsum = zero;
  distsq = zero;
  for (j = 1; j <= n; j++) {
    sum = zero;
    for (k = 1; k <= npt; k++) {
      sum = sum + BMAT(k, j) * W(k);
    }
    jp = j + npt;
    bsum = bsum + sum * W(jp);
    for (ip = npt + 1; ip <= ndim; ip++) {
      sum = sum + BMAT(ip, j) * W(ip);
    }
    bsum = bsum + sum * W(jp);
    VLAG(jp) = sum;
    distsq = distsq + XPT(knew, j) * XPT(knew, j);
  }
  beta = half * distsq * distsq + beta - bsum;
  VLAG(kopt) = VLAG(kopt) + one;
  //
  //     KOLD is set to the index of the provisional interpolation point that
  //     is going to be deleted to make way for the KNEW-th original
  //     interpolation point. The choice of KOLD is governed by the avoidance
  //     of a small value of the denominator in the updating calculation of
  //     UPDATE.
  //
  denom = zero;
  vlmxsq = zero;
  for (k = 1; k <= npt; k++) {
    if (PTSID(k) != zero) {
      hdiag = zero;
      for (j = 1; j <= nptm; j++) {
        hdiag = hdiag + ZMAT(k, j) * ZMAT(k, j);
      }
      den = beta * hdiag + VLAG(k) * VLAG(k);
      if (den > denom) {
        kold = k;
        denom = den;
      }
    }
    vlmxsq = std::max(vlmxsq, VLAG(k) * VLAG(k));
  }
  if (denom <= 1.0e-2 * vlmxsq) {
    W(ndim + knew) = -W(ndim + knew) - winc;
    goto L120;
  }
  goto L80;
  //
  //     When label 260 is reached, all the final positions of the
  //     interpolation points have been chosen although any changes have not
  //     been included yet in XPT. Also the final BMAT and ZMAT matrices are
  //     complete, but, apart from the shift of XBASE, the updating of the
  //     quadratic model remains to be done. The following cycle through the
  //     new interpolation points begins by putting the new point in
  //     XPT(KPT,.) and by setting PQ(KPT) to zero, except that a RETURN
  //     occurs if MAXFUN prohibits another value of F.
  //
L260:
  for (kpt = 1; kpt <= npt; kpt++) {
    if (PTSID(kpt) == zero) continue;  // GOTO 340
    if (nf >= maxfun) {
      nf = -1;
      return;  // GOTO 350
    }
    ih = 0;
    for (j = 1; j <= n; j++) {
      W(j) = XPT(kpt, j);
      XPT(kpt, j) = zero;
      temp = PQ(kpt) * W(j);
      for (i = 1; i <= j; i++) {
        ih = ih + 1;
        HQ(ih) = HQ(ih) + temp * W(i);
      }
    }
    PQ(kpt) = zero;
    ip = static_cast<int>(PTSID(kpt));
    iq = static_cast<int>(static_cast<double>(np) * PTSID(kpt) -
                          static_cast<double>(ip * np));
    if (ip > 0) {
      xp = PTSAUX(1, ip);
      XPT(kpt, ip) = xp;
    }
    if (iq > 0) {
      xq = PTSAUX(1, iq);
      if (ip == 0) xq = PTSAUX(2, iq);
      XPT(kpt, iq) = xq;
    }
    //
    //     Set VQUAD to the value of the current model at the new point.
    //
    vquad = fbase;
    if (ip > 0) {
      ihp = (ip + ip * ip) / 2;
      vquad = vquad + xp * (GOPT(ip) + half * xp * HQ(ihp));
    }
    if (iq > 0) {
      ihq = (iq + iq * iq) / 2;
      vquad = vquad + xq * (GOPT(iq) + half * xq * HQ(ihq));
      if (ip > 0) {
        iw = std::max(ihp, ihq) - std::abs(ip - iq);
        vquad = vquad + xp * xq * HQ(iw);
      }
    }
    for (k = 1; k <= npt; k++) {
      temp = zero;
      if (ip > 0) temp = temp + xp * XPT(k, ip);
      if (iq > 0) temp = temp + xq * XPT(k, iq);
      vquad = vquad + half * PQ(k) * temp * temp;
    }
    //
    //     Calculate F at the new interpolation point, and set DIFF to the
    //     factor that is going to multiply the KPT-th Lagrange function when
    //     the model is updated to provide interpolation to the new function
    //     value.
    //
    for (i = 1; i <= n; i++) {
      W(i) = std::min(std::max(XL(i), XBASE(i) + XPT(kpt, i)), XU(i));
      if (XPT(kpt, i) == SL(i)) W(i) = XL(i);
      if (XPT(kpt, i) == SU(i)) W(i) = XU(i);
    }
    nf = nf + 1;
    // minqa: F = CALFUN (N,X,IPRINT) with X undefined; Powell: CALFUN(N,W,F)
    f = calfun(w);
    FVAL(kpt) = f;
    if (f < FVAL(kopt)) kopt = kpt;
    diff = f - vquad;
    //
    //     Update the quadratic model. The RETURN from the subroutine occurs
    //     when all the new interpolation points are included in the model.
    //
    for (i = 1; i <= n; i++) {
      GOPT(i) = GOPT(i) + diff * BMAT(kpt, i);
    }
    for (k = 1; k <= npt; k++) {
      sum = zero;
      for (j = 1; j <= nptm; j++) {
        sum = sum + ZMAT(k, j) * ZMAT(kpt, j);
      }
      temp = diff * sum;
      if (PTSID(k) == zero) {
        PQ(k) = PQ(k) + temp;
      } else {
        ip = static_cast<int>(PTSID(k));
        iq = static_cast<int>(static_cast<double>(np) * PTSID(k) -
                              static_cast<double>(ip * np));
        ihq = (iq * iq + iq) / 2;
        if (ip == 0) {
          HQ(ihq) = HQ(ihq) + temp * (PTSAUX(2, iq) * PTSAUX(2, iq));
        } else {
          ihp = (ip * ip + ip) / 2;
          HQ(ihp) = HQ(ihp) + temp * (PTSAUX(1, ip) * PTSAUX(1, ip));
          if (iq > 0) {
            HQ(ihq) = HQ(ihq) + temp * (PTSAUX(1, iq) * PTSAUX(1, iq));
            iw = std::max(ihp, ihq) - std::abs(iq - ip);
            HQ(iw) = HQ(iw) + temp * PTSAUX(1, ip) * PTSAUX(1, iq);
          }
        }
      }
    }
    PTSID(kpt) = zero;
  }
}

// ---------------------------------------------------------------------------
// bobyqb.f: SUBROUTINE BOBYQB
// ---------------------------------------------------------------------------
void bobyqb(int n, int npt, double *x, const double *xl, const double *xu,
            double rhobeg, double rhoend, int maxfun, double *xbase,
            double *xpt, double *fval, double *xopt, double *gopt, double *hq,
            double *pq, double *bmat, double *zmat, int ndim, double *sl,
            double *su, double *xnew, double *xalt, double *d, double *vlag,
            double *w, int &ierr, MinqaCalfun &calfun) {
  double half, one, ten, tenth, two, zero, xoptsq, fsave, rho, delta, diffa,
      diffb, diffc = 0.0, temp, dsq = 0.0, crvmin = 0.0, dnorm = 0.0, distsq,
      errbig, frhosq, bdtol, bdtest, curv, fracsq, sumpq, sum, sumz, sumw,
      adelt = 0.0, alpha = 0.0, cauchy = 0.0, suma, sumb, beta = 0.0, bsum,
      dx, denom = 0.0, delsq, scaden, biglsq, hdiag, den, f = 0.0, fopt = 0.0,
      vquad = 0.0, diff = 0.0, ratio = 0.0, densav, pqold, gqsq, gisq, dist;
  int np, nptm, nh, nf = 0, kopt = 0, i, kbase, nresc, ntrits, itest, nfsav,
      ih, j, k, jj, ip, jp, knew = 0, ksav;

  half = 0.5;
  one = 1.0;
  ten = 10.0;
  tenth = 0.1;
  two = 2.0;
  zero = 0.0;
  np = n + 1;
  nptm = npt - np;
  nh = (n * np) / 2;
  //
  //     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ,
  //     PQ, BMAT and ZMAT for the first iteration, with the corresponding
  //     values of of NF and KOPT, which are the number of calls of CALFUN so
  //     far and the index of the interpolation point at the trust region
  //     centre. Then the initial XOPT is set too. The branch to label 720
  //     occurs if MAXFUN is less than NPT. GOPT will be updated if KOPT is
  //     different from KBASE.
  //
  prelim(n, npt, x, xl, xu, rhobeg, maxfun, xbase, xpt, fval, gopt, hq, pq,
         bmat, zmat, ndim, sl, su, nf, kopt, calfun);
  xoptsq = zero;
  for (i = 1; i <= n; i++) {
    XOPT(i) = XPT(kopt, i);
    xoptsq = xoptsq + XOPT(i) * XOPT(i);
  }
  fsave = FVAL(1);
  if (nf < npt) {
    ierr = 390;
    goto L720;
  }
  kbase = 1;
  //
  //     Complete the settings that are required for the iterative procedure.
  //
  rho = rhobeg;
  delta = rho;
  nresc = nf;
  ntrits = 0;
  diffa = zero;
  diffb = zero;
  itest = 0;
  nfsav = nf;
  //
  //     Update GOPT if necessary before the first iteration and after each
  //     call of RESCUE that makes a call of CALFUN.
  //
L20:
  if (kopt != kbase) {
    ih = 0;
    for (j = 1; j <= n; j++) {
      for (i = 1; i <= j; i++) {
        ih = ih + 1;
        if (i < j) GOPT(j) = GOPT(j) + HQ(ih) * XOPT(i);
        GOPT(i) = GOPT(i) + HQ(ih) * XOPT(j);
      }
    }
    if (nf > npt) {
      for (k = 1; k <= npt; k++) {
        temp = zero;
        for (j = 1; j <= n; j++) {
          temp = temp + XPT(k, j) * XOPT(j);
        }
        temp = PQ(k) * temp;
        for (i = 1; i <= n; i++) {
          GOPT(i) = GOPT(i) + temp * XPT(k, i);
        }
      }
    }
  }
  //
  //     Generate the next point in the trust region that provides a small
  //     value of the quadratic model subject to the constraints on the
  //     variables. The integer NTRITS is set to the number "trust region"
  //     iterations that have occurred since the last "alternative" iteration.
  //     If the length of XNEW-XOPT is less than HALF*RHO, however, then there
  //     is a branch to label 650 or 680 with NTRITS=-1, instead of
  //     calculating F at XNEW.
  //
L60:
  trsbox(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, xnew, d, w,
         &W(np), &W(np + n), &W(np + 2 * n), &W(np + 3 * n), dsq, crvmin);
  dnorm = std::min(delta, std::sqrt(dsq));
  if (dnorm < half * rho) {
    ntrits = -1;
    distsq = (ten * rho) * (ten * rho);
    if (nf <= nfsav + 2) goto L650;
    //
    //     The following choice between labels 650 and 680 depends on whether
    //     or not our work with the current RHO seems to be complete. Either
    //     RHO is decreased or termination occurs if the errors in the
    //     quadratic model at the last three interpolation points compare
    //     favourably with predictions of likely improvements to the model
    //     within distance HALF*RHO of XOPT.
    //
    errbig = std::max(std::max(diffa, diffb), diffc);
    frhosq = 0.125 * rho * rho;
    if (crvmin > zero && errbig > frhosq * crvmin) goto L650;
    bdtol = errbig / rho;
    for (j = 1; j <= n; j++) {
      bdtest = bdtol;
      if (XNEW(j) == SL(j)) bdtest = W(j);
      if (XNEW(j) == SU(j)) bdtest = -W(j);
      if (bdtest < bdtol) {
        curv = HQ((j + j * j) / 2);
        for (k = 1; k <= npt; k++) {
          curv = curv + PQ(k) * (XPT(k, j) * XPT(k, j));
        }
        bdtest = bdtest + half * curv * rho;
        if (bdtest < bdtol) goto L650;
      }
    }
    goto L680;
  }
  ntrits = ntrits + 1;
  //
  //     Severe cancellation is likely to occur if XOPT is too far from XBASE.
  //     If the following test holds, then XBASE is shifted so that XOPT
  //     becomes zero. The appropriate changes are made to BMAT and to the
  //     second derivatives of the current model, beginning with the changes
  //     to BMAT that do not depend on ZMAT. VLAG is used temporarily for
  //     working space.
  //
L90:
  if (dsq <= 1.0e-3 * xoptsq) {
    fracsq = 0.25 * xoptsq;
    sumpq = zero;
    for (k = 1; k <= npt; k++) {
      sumpq = sumpq + PQ(k);
      sum = -half * xoptsq;
      for (i = 1; i <= n; i++) {
        sum = sum + XPT(k, i) * XOPT(i);
      }
      W(npt + k) = sum;
      temp = fracsq - half * sum;
      for (i = 1; i <= n; i++) {
        W(i) = BMAT(k, i);
        VLAG(i) = sum * XPT(k, i) + temp * XOPT(i);
        ip = npt + i;
        for (j = 1; j <= i; j++) {
          BMAT(ip, j) = BMAT(ip, j) + W(i) * VLAG(j) + VLAG(i) * W(j);
        }
      }
    }
    //
    //     Then the revisions of BMAT that depend on ZMAT are calculated.
    //
    for (jj = 1; jj <= nptm; jj++) {
      sumz = zero;
      sumw = zero;
      for (k = 1; k <= npt; k++) {
        sumz = sumz + ZMAT(k, jj);
        VLAG(k) = W(npt + k) * ZMAT(k, jj);
        sumw = sumw + VLAG(k);
      }
      for (j = 1; j <= n; j++) {
        sum = (fracsq * sumz - half * sumw) * XOPT(j);
        for (k = 1; k <= npt; k++) {
          sum = sum + VLAG(k) * XPT(k, j);
        }
        W(j) = sum;
        for (k = 1; k <= npt; k++) {
          BMAT(k, j) = BMAT(k, j) + sum * ZMAT(k, jj);
        }
      }
      for (i = 1; i <= n; i++) {
        ip = i + npt;
        temp = W(i);
        for (j = 1; j <= i; j++) {
          BMAT(ip, j) = BMAT(ip, j) + temp * W(j);
        }
      }
    }
    //
    //     The following instructions complete the shift, including the
    //     changes to the second derivative parameters of the quadratic model.
    //
    ih = 0;
    for (j = 1; j <= n; j++) {
      W(j) = -half * sumpq * XOPT(j);
      for (k = 1; k <= npt; k++) {
        W(j) = W(j) + PQ(k) * XPT(k, j);
        XPT(k, j) = XPT(k, j) - XOPT(j);
      }
      for (i = 1; i <= j; i++) {
        ih = ih + 1;
        HQ(ih) = HQ(ih) + W(i) * XOPT(j) + XOPT(i) * W(j);
        BMAT(npt + i, j) = BMAT(npt + j, i);
      }
    }
    for (i = 1; i <= n; i++) {
      XBASE(i) = XBASE(i) + XOPT(i);
      XNEW(i) = XNEW(i) - XOPT(i);
      SL(i) = SL(i) - XOPT(i);
      SU(i) = SU(i) - XOPT(i);
      XOPT(i) = zero;
    }
    xoptsq = zero;
  }
  if (ntrits == 0) goto L210;
  goto L230;
  //
  //     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
  //     more expensive than the previous shift, because new matrices BMAT and
  //     ZMAT are generated from scratch, which may include the replacement of
  //     interpolation points whose positions seem to be causing near linear
  //     dependence in the interpolation conditions. Therefore RESCUE is
  //     called only if rounding errors have reduced by at least a factor of
  //     two the denominator of the formula for updating the H matrix. It
  //     provides a useful safeguard, but is not invoked in most applications
  //     of BOBYQA.
  //
L190:
  nfsav = nf;
  kbase = kopt;
  rescue(n, npt, xl, xu, maxfun, xbase, xpt, fval, xopt, gopt, hq, pq, bmat,
         zmat, ndim, sl, su, nf, delta, kopt, vlag, w, &W(n + np),
         &W(ndim + np), calfun);
  //
  //     XOPT is updated now in case the branch below to label 720 is taken.
  //     Any updating of GOPT occurs after the branch below to label 20, which
  //     leads to a trust region iteration as does the branch to label 60.
  //
  xoptsq = zero;
  if (kopt != kbase) {
    for (i = 1; i <= n; i++) {
      XOPT(i) = XPT(kopt, i);
      xoptsq = xoptsq + XOPT(i) * XOPT(i);
    }
  }
  if (nf < 0) {
    nf = maxfun;
    ierr = 390;
    goto L720;
  }
  nresc = nf;
  if (nfsav < nf) {
    nfsav = nf;
    goto L20;
  }
  if (ntrits > 0) goto L60;
  //
  //     Pick two alternative vectors of variables, relative to XBASE, that
  //     are suitable as new positions of the KNEW-th interpolation point.
  //     Firstly, XNEW is set to the point on a line through XOPT and another
  //     interpolation point that minimizes the predicted value of the next
  //     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
  //     and SU bounds. Secondly, XALT is set to the best feasible point on
  //     a constrained version of the Cauchy step of the KNEW-th Lagrange
  //     function, the corresponding value of the square of this function
  //     being returned in CAUCHY. The choice between these alternatives is
  //     going to be made when the denominator is calculated.
  //
L210:
  altmov(n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, knew, adelt, xnew,
         xalt, alpha, cauchy, w, &W(np), &W(ndim + 1));
  for (i = 1; i <= n; i++) {
    D(i) = XNEW(i) - XOPT(i);
  }
  //
  //     Calculate VLAG and BETA for the current choice of D. The scalar
  //     product of D with XPT(K,.) is going to be held in W(NPT+K) for
  //     use when VQUAD is calculated.
  //
L230:
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
    W(npt + k) = suma;
  }
  beta = zero;
  for (jj = 1; jj <= nptm; jj++) {
    sum = zero;
    for (k = 1; k <= npt; k++) {
      sum = sum + ZMAT(k, jj) * W(k);
    }
    beta = beta - sum * sum;
    for (k = 1; k <= npt; k++) {
      VLAG(k) = VLAG(k) + sum * ZMAT(k, jj);
    }
  }
  dsq = zero;
  bsum = zero;
  dx = zero;
  for (j = 1; j <= n; j++) {
    dsq = dsq + D(j) * D(j);
    sum = zero;
    for (k = 1; k <= npt; k++) {
      sum = sum + W(k) * BMAT(k, j);
    }
    bsum = bsum + sum * D(j);
    jp = npt + j;
    for (i = 1; i <= n; i++) {
      sum = sum + BMAT(jp, i) * D(i);
    }
    VLAG(jp) = sum;
    bsum = bsum + sum * D(j);
    dx = dx + D(j) * XOPT(j);
  }
  beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
  VLAG(kopt) = VLAG(kopt) + one;
  //
  //     If NTRITS is zero, the denominator may be increased by replacing
  //     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
  //     rounding errors have damaged the chosen denominator.
  //
  if (ntrits == 0) {
    denom = VLAG(knew) * VLAG(knew) + alpha * beta;
    if (denom < cauchy && cauchy > zero) {
      for (i = 1; i <= n; i++) {
        XNEW(i) = XALT(i);
        D(i) = XNEW(i) - XOPT(i);
      }
      cauchy = zero;
      goto L230;
    }
    if (denom <= half * (VLAG(knew) * VLAG(knew))) {
      if (nf > nresc) goto L190;
      ierr = 320;
      goto L720;
    }
    //
    //     Alternatively, if NTRITS is positive, then set KNEW to the index of
    //     the next interpolation point to be deleted to make room for a trust
    //     region step. Again RESCUE may be called if rounding errors have
    //     damaged the chosen denominator, which is the reason for attempting
    //     to select KNEW before calculating the next value of the objective
    //     function.
    //
  } else {
    delsq = delta * delta;
    scaden = zero;
    biglsq = zero;
    knew = 0;
    for (k = 1; k <= npt; k++) {
      if (k == kopt) continue;  // GOTO 350
      hdiag = zero;
      for (jj = 1; jj <= nptm; jj++) {
        hdiag = hdiag + ZMAT(k, jj) * ZMAT(k, jj);
      }
      den = beta * hdiag + VLAG(k) * VLAG(k);
      distsq = zero;
      for (j = 1; j <= n; j++) {
        distsq = distsq + (XPT(k, j) - XOPT(j)) * (XPT(k, j) - XOPT(j));
      }
      temp = std::max(one, (distsq / delsq) * (distsq / delsq));
      if (temp * den > scaden) {
        scaden = temp * den;
        knew = k;
        denom = den;
      }
      biglsq = std::max(biglsq, temp * (VLAG(k) * VLAG(k)));
    }
    if (scaden <= half * biglsq) {
      if (nf > nresc) goto L190;
      ierr = 320;
      goto L720;
    }
  }
  //
  //     Put the variables for the next calculation of the objective function
  //       in XNEW, with any adjustments for the bounds.
  //
  //
  //     Calculate the value of the objective function at XBASE+XNEW, unless
  //       the limit on the number of calculations of F has been reached.
  //
L360:
  for (i = 1; i <= n; i++) {
    X(i) = std::min(std::max(XL(i), XBASE(i) + XNEW(i)), XU(i));
    if (XNEW(i) == SL(i)) X(i) = XL(i);
    if (XNEW(i) == SU(i)) X(i) = XU(i);
  }
  if (nf >= maxfun) {
    ierr = 390;
    goto L720;
  }
  nf = nf + 1;
  f = calfun(x);
  if (ntrits == -1) {
    fsave = f;
    goto L720;
  }
  //
  //     Use the quadratic model to predict the change in F due to the step D,
  //       and set DIFF to the error of this prediction.
  //
  fopt = FVAL(kopt);
  vquad = zero;
  ih = 0;
  for (j = 1; j <= n; j++) {
    vquad = vquad + D(j) * GOPT(j);
    for (i = 1; i <= j; i++) {
      ih = ih + 1;
      temp = D(i) * D(j);
      if (i == j) temp = half * temp;
      vquad = vquad + HQ(ih) * temp;
    }
  }
  for (k = 1; k <= npt; k++) {
    vquad = vquad + half * PQ(k) * (W(npt + k) * W(npt + k));
  }
  diff = f - fopt - vquad;
  diffc = diffb;
  diffb = diffa;
  diffa = std::fabs(diff);
  if (dnorm > rho) nfsav = nf;
  //
  //     Pick the next value of DELTA after a trust region step.
  //
  if (ntrits > 0) {
    if (vquad >= zero) {
      ierr = 430;
      goto L720;
    }
    ratio = (f - fopt) / vquad;
    if (ratio <= tenth) {
      delta = std::min(half * delta, dnorm);
    } else if (ratio <= 0.7) {
      delta = std::max(half * delta, dnorm);
    } else {
      delta = std::max(half * delta, dnorm + dnorm);
    }
    if (delta <= 1.5 * rho) delta = rho;
    //
    //     Recalculate KNEW and DENOM if the new F is less than FOPT.
    //
    if (f < fopt) {
      ksav = knew;
      densav = denom;
      delsq = delta * delta;
      scaden = zero;
      biglsq = zero;
      knew = 0;
      for (k = 1; k <= npt; k++) {
        hdiag = zero;
        for (jj = 1; jj <= nptm; jj++) {
          hdiag = hdiag + ZMAT(k, jj) * ZMAT(k, jj);
        }
        den = beta * hdiag + VLAG(k) * VLAG(k);
        distsq = zero;
        for (j = 1; j <= n; j++) {
          distsq = distsq + (XPT(k, j) - XNEW(j)) * (XPT(k, j) - XNEW(j));
        }
        temp = std::max(one, (distsq / delsq) * (distsq / delsq));
        if (temp * den > scaden) {
          scaden = temp * den;
          knew = k;
          denom = den;
        }
        biglsq = std::max(biglsq, temp * (VLAG(k) * VLAG(k)));
      }
      if (scaden <= half * biglsq) {
        knew = ksav;
        denom = densav;
      }
    }
  }
  //
  //     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
  //     moved. Also update the second derivative terms of the model.
  //
  updatebobyqa(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w);
  ih = 0;
  pqold = PQ(knew);
  PQ(knew) = zero;
  for (i = 1; i <= n; i++) {
    temp = pqold * XPT(knew, i);
    for (j = 1; j <= i; j++) {
      ih = ih + 1;
      HQ(ih) = HQ(ih) + temp * XPT(knew, j);
    }
  }
  for (jj = 1; jj <= nptm; jj++) {
    temp = diff * ZMAT(knew, jj);
    for (k = 1; k <= npt; k++) {
      PQ(k) = PQ(k) + temp * ZMAT(k, jj);
    }
  }
  //
  //     Include the new interpolation point, and make the changes to GOPT at
  //     the old XOPT that are caused by the updating of the quadratic model.
  //
  FVAL(knew) = f;
  for (i = 1; i <= n; i++) {
    XPT(knew, i) = XNEW(i);
    W(i) = BMAT(knew, i);
  }
  for (k = 1; k <= npt; k++) {
    suma = zero;
    for (jj = 1; jj <= nptm; jj++) {
      suma = suma + ZMAT(knew, jj) * ZMAT(k, jj);
    }
    sumb = zero;
    for (j = 1; j <= n; j++) {
      sumb = sumb + XPT(k, j) * XOPT(j);
    }
    temp = suma * sumb;
    for (i = 1; i <= n; i++) {
      W(i) = W(i) + temp * XPT(k, i);
    }
  }
  for (i = 1; i <= n; i++) {
    GOPT(i) = GOPT(i) + diff * W(i);
  }
  //
  //     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
  //
  if (f < fopt) {
    kopt = knew;
    xoptsq = zero;
    ih = 0;
    for (j = 1; j <= n; j++) {
      XOPT(j) = XNEW(j);
      xoptsq = xoptsq + XOPT(j) * XOPT(j);
      for (i = 1; i <= j; i++) {
        ih = ih + 1;
        if (i < j) GOPT(j) = GOPT(j) + HQ(ih) * D(i);
        GOPT(i) = GOPT(i) + HQ(ih) * D(j);
      }
    }
    for (k = 1; k <= npt; k++) {
      temp = zero;
      for (j = 1; j <= n; j++) {
        temp = temp + XPT(k, j) * D(j);
      }
      temp = PQ(k) * temp;
      for (i = 1; i <= n; i++) {
        GOPT(i) = GOPT(i) + temp * XPT(k, i);
      }
    }
  }
  //
  //     Calculate the parameters of the least Frobenius norm interpolant to
  //     the current data, the gradient of this interpolant at XOPT being put
  //     into VLAG(NPT+I), I=1,2,...,N.
  //
  if (ntrits > 0) {
    for (k = 1; k <= npt; k++) {
      VLAG(k) = FVAL(k) - FVAL(kopt);
      W(k) = zero;
    }
    for (j = 1; j <= nptm; j++) {
      sum = zero;
      for (k = 1; k <= npt; k++) {
        sum = sum + ZMAT(k, j) * VLAG(k);
      }
      for (k = 1; k <= npt; k++) {
        W(k) = W(k) + sum * ZMAT(k, j);
      }
    }
    for (k = 1; k <= npt; k++) {
      sum = zero;
      for (j = 1; j <= n; j++) {
        sum = sum + XPT(k, j) * XOPT(j);
      }
      W(k + npt) = W(k);
      W(k) = sum * W(k);
    }
    gqsq = zero;
    gisq = zero;
    for (i = 1; i <= n; i++) {
      sum = zero;
      for (k = 1; k <= npt; k++) {
        sum = sum + BMAT(k, i) * VLAG(k) + XPT(k, i) * W(k);
      }
      if (XOPT(i) == SL(i)) {
        gqsq = gqsq + std::min(zero, GOPT(i)) * std::min(zero, GOPT(i));
        gisq = gisq + std::min(zero, sum) * std::min(zero, sum);
      } else if (XOPT(i) == SU(i)) {
        gqsq = gqsq + std::max(zero, GOPT(i)) * std::max(zero, GOPT(i));
        gisq = gisq + std::max(zero, sum) * std::max(zero, sum);
      } else {
        gqsq = gqsq + GOPT(i) * GOPT(i);
        gisq = gisq + sum * sum;
      }
      VLAG(npt + i) = sum;
    }
    //
    //     Test whether to replace the new quadratic model by the least
    //     Frobenius norm interpolant, making the replacement if the test is
    //     satisfied.
    //
    itest = itest + 1;
    if (gqsq < ten * gisq) itest = 0;
    if (itest >= 3) {
      for (i = 1; i <= std::max(npt, nh); i++) {
        if (i <= n) GOPT(i) = VLAG(npt + i);
        if (i <= npt) PQ(i) = W(npt + i);
        if (i <= nh) HQ(i) = zero;
        itest = 0;
      }
    }
  }
  //
  //     If a trust region step has provided a sufficient decrease in F, then
  //     branch for another trust region calculation. The case NTRITS=0
  //     occurs when the new interpolation point was reached by an alternative
  //     step.
  //
  if (ntrits == 0) goto L60;
  if (f <= fopt + tenth * vquad) goto L60;
  //
  //     Alternatively, find out if the interpolation points are close enough
  //       to the best point so far.
  //
  distsq = std::max((two * delta) * (two * delta), (ten * rho) * (ten * rho));
L650:
  knew = 0;
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
  //
  //     If KNEW is positive, then ALTMOV finds alternative new positions for
  //     the KNEW-th interpolation point within distance ADELT of XOPT. It is
  //     reached via label 90. Otherwise, there is a branch to label 60 for
  //     another trust region iteration, unless the calculations with the
  //     current RHO are complete.
  //
  if (knew > 0) {
    dist = std::sqrt(distsq);
    if (ntrits == -1) {
      delta = std::min(tenth * delta, half * dist);
      if (delta <= 1.5 * rho) delta = rho;
    }
    ntrits = 0;
    adelt = std::max(std::min(tenth * dist, delta), rho);
    dsq = adelt * adelt;
    goto L90;
  }
  if (ntrits == -1) goto L680;
  if (ratio > zero) goto L60;
  if (std::max(delta, dnorm) > rho) goto L60;
  //
  //     The calculations with the current value of RHO are complete. Pick the
  //       next values of RHO and DELTA.
  //
L680:
  if (rho > rhoend) {
    delta = half * rho;
    ratio = rho / rhoend;
    if (ratio <= 16.0) {
      rho = rhoend;
    } else if (ratio <= 250.0) {
      rho = std::sqrt(ratio) * rhoend;
    } else {
      rho = tenth * rho;
    }
    delta = std::max(delta, rho);
    calfun.minqit(rho, nf, FVAL(kopt), xbase, xopt);
    ntrits = 0;
    nfsav = nf;
    goto L60;
  }
  //
  //     Return from the calculation, after another Newton-Raphson step, if
  //       it is too short to have been tried before.
  //
  if (ntrits == -1) goto L360;
L720:
  if (FVAL(kopt) <= fsave) {
    for (i = 1; i <= n; i++) {
      X(i) = std::min(std::max(XL(i), XBASE(i) + XOPT(i)), XU(i));
      if (XOPT(i) == SL(i)) X(i) = XL(i);
      if (XOPT(i) == SU(i)) X(i) = XU(i);
    }
    f = FVAL(kopt);
  }
  calfun.minqir(f, nf, x);
}

// ---------------------------------------------------------------------------
// bobyqa.f: SUBROUTINE BOBYQA
// ---------------------------------------------------------------------------
void bobyqa(int n, int npt, double *x, const double *xl, const double *xu,
            double rhobeg, double rhoend, int maxfun, double *w, int &ierr,
            MinqaCalfun &calfun) {
  double zero, temp;
  int np, ndim, ixb, ixp, ifv, ixo, igo, ihq, ipq, ibmat, izmat, isl, isu,
      ixn, ixa, id, ivl, iw, j, jsl, jsu;

  np = n + 1;
  if (npt < n + 2 || npt > ((n + 2) * np) / 2) {
    ierr = 10;
    return;  // GO TO 40
  }
  //
  //     Partition the working space array, so that different parts of it can
  //     be treated separately during the calculation of BOBYQB. The partition
  //     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
  //     space that is taken by the last array in the argument list of BOBYQB.
  //
  ndim = npt + n;
  ixb = 1;
  ixp = ixb + n;
  ifv = ixp + n * npt;
  ixo = ifv + npt;
  igo = ixo + n;
  ihq = igo + n;
  ipq = ihq + (n * np) / 2;
  ibmat = ipq + npt;
  izmat = ibmat + ndim * n;
  isl = izmat + npt * (npt - np);
  isu = isl + n;
  ixn = isu + n;
  ixa = ixn + n;
  id = ixa + n;
  ivl = id + n;
  iw = ivl + ndim;
  ierr = 0;
  //
  //     Return if there is insufficient space between the bounds. Modify the
  //     initial X if necessary in order to avoid conflicts between the bounds
  //     and the construction of the first quadratic model. The lower and
  //     upper bounds on moves from the updated X are set now, in the ISL and
  //     ISU partitions of W, in order to provide useful and exact information
  //     about components of X that become within distance RHOBEG from their
  //     bounds.
  //
  zero = 0.0;
  for (j = 1; j <= n; j++) {
    temp = XU(j) - XL(j);
    if (temp < rhobeg + rhobeg) {
      ierr = 20;
      return;  // GOTO 40
    }
    jsl = isl + j - 1;
    jsu = jsl + n;
    W(jsl) = XL(j) - X(j);
    W(jsu) = XU(j) - X(j);
    if (W(jsl) >= -rhobeg) {
      if (W(jsl) >= zero) {
        X(j) = XL(j);
        W(jsl) = zero;
        W(jsu) = temp;
      } else {
        X(j) = XL(j) + rhobeg;
        W(jsl) = -rhobeg;
        W(jsu) = std::max(XU(j) - X(j), rhobeg);
      }
    } else if (W(jsu) <= rhobeg) {
      if (W(jsu) <= zero) {
        X(j) = XU(j);
        W(jsl) = -temp;
        W(jsu) = zero;
      } else {
        X(j) = XU(j) - rhobeg;
        W(jsl) = std::min(XL(j) - X(j), -rhobeg);
        W(jsu) = rhobeg;
      }
    }
  }
  //
  //     Make the call of BOBYQB.
  //
  bobyqb(n, npt, x, xl, xu, rhobeg, rhoend, maxfun, &W(ixb), &W(ixp),
         &W(ifv), &W(ixo), &W(igo), &W(ihq), &W(ipq), &W(ibmat), &W(izmat),
         ndim, &W(isl), &W(isu), &W(ixn), &W(ixa), &W(id), &W(ivl), &W(iw),
         ierr, calfun);
}

}  // namespace

#undef X
#undef XL
#undef XU
#undef W
#undef XBASE
#undef XPT
#undef FVAL
#undef XOPT
#undef GOPT
#undef HQ
#undef PQ
#undef BMAT
#undef ZMAT
#undef SL
#undef SU
#undef XNEW
#undef XALT
#undef D
#undef VLAG
#undef GNEW
#undef XBDI
#undef S
#undef HS
#undef HRED
#undef GLAG
#undef HCOL
#undef PTSAUX
#undef PTSID

// Public entry point (declared in minqa_types.h). Mirrors minqa's
// bobyqa_cpp(): workspace sized as (NPT+5)*(NPT+N)+3*N*(N+5)/2, a copy of
// par is optimized in place, then par/fval are finalized by minqa_run().
// A NULL lower (upper) means -Inf (+Inf) for every component; NULL opts
// means minqa_options_default(n, par).
int bobyqa_solve_c(int n, const double *par, const double *lower,
                   const double *upper, minqa_c_objfun_t objfun,
                   void *userdata, const minqa_options_t *opts,
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
  std::vector<double> x, xl, xu, w;
  try {
    x.assign(par, par + n);
    xl.assign(static_cast<size_t>(n), -HUGE_VAL);
    xu.assign(static_cast<size_t>(n), HUGE_VAL);
    if (lower != nullptr) xl.assign(lower, lower + n);
    if (upper != nullptr) xu.assign(upper, upper + n);
    // minqa.cpp: vector<double> w((np + 5) * (np + n) + (3 * n * (n + 5))/2)
    long long lw = (static_cast<long long>(opts->npt) + 5) *
                       (static_cast<long long>(opts->npt) + n) +
                   (3LL * n * (n + 5)) / 2;
    if (lw < 1) lw = 1;
    w.assign(static_cast<size_t>(lw), 0.0);
  } catch (const std::bad_alloc &) {
    minqa_result_zero(result);
    result->n = n;
    result->ierr = MINQA_ERR_NOMEM;
    return MINQA_ERR_NOMEM;
  }
  return minqa_run(n, x.data(), calfun, result, [&]() -> int {
    int ierr = 0;
    bobyqa(n, opts->npt, x.data(), xl.data(), xu.data(), opts->rhobeg,
           opts->rhoend, opts->maxfun, w.data(), ierr, calfun);
    return ierr;
  });
}
