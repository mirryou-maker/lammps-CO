// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   OpenMP-threaded variant of pair_lj_cut_dipole_long.
   Follows the ThrOMP pattern established in pair_lj_cut_dipole_cut_omp.
   Torque on outer atom i accumulated in t1tmp/t2tmp/t3tmp (thread-safe).
   Torque on inner atom j uses per-thread thr->get_torque() buffer.
------------------------------------------------------------------------- */

#include "pair_lj_cut_dipole_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

namespace { using dbl3_t = struct { double x,y,z; }; }

/* ---------------------------------------------------------------------- */

PairLJCutDipoleLongOMP::PairLJCutDipoleLongOMP(LAMMPS *lmp) :
  PairLJCutDipoleLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJCutDipoleLongOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else                    eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else                    eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else                    eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutDipoleLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * const * const torque = thr->get_torque();
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj   = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int ** const firstneigh             = list->firstneigh;

  const double pre1_factor = 2.0 * g_ewald / MY_PIS;
  const double pre2_factor = 4.0 * g_ewald*g_ewald*g_ewald / MY_PIS;
  const double pre3_factor = 8.0 * g_ewald*g_ewald*g_ewald*g_ewald*g_ewald / MY_PIS;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    const double qtmp  = q[i];
    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const int itype    = type[i];
    const double * _noalias const cutsqi    = cutsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const lj1i      = lj1[itype];
    const double * _noalias const lj2i      = lj2[itype];
    const double * _noalias const lj3i      = lj3[itype];
    const double * _noalias const lj4i      = lj4[itype];
    const double * _noalias const offseti   = offset[itype];
    const double mui0 = atom->mu[i][0];
    const double mui1 = atom->mu[i][1];
    const double mui2 = atom->mu[i][2];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;
    double t1tmp = 0.0, t2tmp = 0.0, t3tmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_lj   = special_lj[sbmask(j)];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq  = delx*delx + dely*dely + delz*delz;
      const int jtype   = type[j];

      if (rsq < cutsqi[jtype]) {
        const double r2inv = 1.0/rsq;
        const double rinv  = sqrt(r2inv);

        double forcecoulx = 0.0, forcecouly = 0.0, forcecoulz = 0.0;
        double tixcoul = 0.0, tiycoul = 0.0, tizcoul = 0.0;
        double tjxcoul = 0.0, tjycoul = 0.0, tjzcoul = 0.0;
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
        double g0 = 0.0, g1 = 0.0, g2 = 0.0;
        double erfc_val = 0.0, expm2_val = 0.0;
        double pdotp = 0.0, pidotr = 0.0, pjdotr = 0.0;
        const double muj0 = atom->mu[j][0];
        const double muj1 = atom->mu[j][1];
        const double muj2 = atom->mu[j][2];

        if (rsq < cut_coulsq) {
          const double r    = sqrt(rsq);
          const double grij = g_ewald * r;
          expm2_val = exp(-grij*grij);
          const double t = 1.0 / (1.0 + EWALD_P*grij);
          erfc_val = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2_val;

          pdotp = mui0*muj0 + mui1*muj1 + mui2*muj2;
          pidotr = mui0*delx + mui1*dely + mui2*delz;
          pjdotr = muj0*delx + muj1*dely + muj2*delz;

          g0 = qtmp*q[j];
          g1 = qtmp*pjdotr - q[j]*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > 0.0) {
            b0 = erfc_val * rinv;
            b1 = (b0 + pre1_factor*expm2_val) * r2inv;
            b2 = (3.0*b1 + pre2_factor*expm2_val) * r2inv;
            const double b3 = (5.0*b2 + pre3_factor*expm2_val) * r2inv;

            const double g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            double fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj0 - q[j]*mui0) +
              b2 * (pjdotr*mui0 + pidotr*muj0);
            double fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj1 - q[j]*mui1) +
              b2 * (pjdotr*mui1 + pidotr*muj1);
            double fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj2 - q[j]*mui2) +
              b2 * (pjdotr*mui2 + pidotr*muj2);

            double zdix = delx * (q[j]*b1 + b2*pjdotr) - b1*muj0;
            double zdiy = dely * (q[j]*b1 + b2*pjdotr) - b1*muj1;
            double zdiz = delz * (q[j]*b1 + b2*pjdotr) - b1*muj2;
            double zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mui0;
            double zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mui1;
            double zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mui2;

            if (factor_coul < 1.0) {
              fdx *= factor_coul; fdy *= factor_coul; fdz *= factor_coul;
              zdix *= factor_coul; zdiy *= factor_coul; zdiz *= factor_coul;
              zdjx *= factor_coul; zdjy *= factor_coul; zdjz *= factor_coul;
            }

            forcecoulx += fdx; forcecouly += fdy; forcecoulz += fdz;
            tixcoul += mui1*zdiz - mui2*zdiy;
            tiycoul += mui2*zdix - mui0*zdiz;
            tizcoul += mui0*zdiy - mui1*zdix;
            tjxcoul += muj1*zdjz - muj2*zdjy;
            tjycoul += muj2*zdjx - muj0*zdjz;
            tjzcoul += muj0*zdjy - muj1*zdjx;
          }

          if (factor_coul < 1.0) {
            d0 = (erfc_val - 1.0) * rinv;
            d1 = (d0 + pre1_factor*expm2_val) * r2inv;
            d2 = (3.0*d1 + pre2_factor*expm2_val) * r2inv;
            d3 = (5.0*d2 + pre3_factor*expm2_val) * r2inv;

            const double g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            double fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj0 - q[j]*mui0) +
              d2 * (pjdotr*mui0 + pidotr*muj0);
            double fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj1 - q[j]*mui1) +
              d2 * (pjdotr*mui1 + pidotr*muj1);
            double faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj2 - q[j]*mui2) +
              d2 * (pjdotr*mui2 + pidotr*muj2);

            double zaix = delx * (q[j]*d1 + d2*pjdotr) - d1*muj0;
            double zaiy = dely * (q[j]*d1 + d2*pjdotr) - d1*muj1;
            double zaiz = delz * (q[j]*d1 + d2*pjdotr) - d1*muj2;
            double zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mui0;
            double zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mui1;
            double zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mui2;

            if (factor_coul > 0.0) {
              const double facm1 = 1.0 - factor_coul;
              fax *= facm1; fay *= facm1; faz *= facm1;
              zaix *= facm1; zaiy *= facm1; zaiz *= facm1;
              zajx *= facm1; zajy *= facm1; zajz *= facm1;
            }

            forcecoulx += fax; forcecouly += fay; forcecoulz += faz;
            tixcoul += mui1*zaiz - mui2*zaiy;
            tiycoul += mui2*zaix - mui0*zaiz;
            tizcoul += mui0*zaiy - mui1*zaix;
            tjxcoul += muj1*zajz - muj2*zajy;
            tjycoul += muj2*zajx - muj0*zajz;
            tjzcoul += muj0*zajy - muj1*zajx;
          }
        }

        double fforce = 0.0;
        double r6inv = 0.0;
        if (rsq < cut_ljsqi[jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          fforce = r6inv * (lj1i[jtype]*r6inv - lj2i[jtype]) * factor_lj * r2inv;
        }

        const double fx = qqrd2e*forcecoulx + delx*fforce;
        const double fy = qqrd2e*forcecouly + dely*fforce;
        const double fz = qqrd2e*forcecoulz + delz*fforce;

        fxtmp += fx;
        fytmp += fy;
        fztmp += fz;
        t1tmp += qqrd2e*tixcoul;
        t2tmp += qqrd2e*tiycoul;
        t3tmp += qqrd2e*tizcoul;

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
          torque[j][0] += qqrd2e*tjxcoul;
          torque[j][1] += qqrd2e*tjycoul;
          torque[j][2] += qqrd2e*tjzcoul;
        }

        double evdwl = 0.0, ecoul = 0.0;
        if (EFLAG) {
          if (rsq < cut_coulsq && factor_coul > 0.0) {
            ecoul = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < 1.0) {
              ecoul *= factor_coul;
              ecoul += (1.0-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          }
          if (rsq < cut_ljsqi[jtype]) {
            evdwl = r6inv*(lj3i[jtype]*r6inv-lj4i[jtype]) - offseti[jtype];
            evdwl *= factor_lj;
          }
        }

        if (EVFLAG) ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR,
                                      evdwl, ecoul, fx, fy, fz, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    torque[i][0] += t1tmp;
    torque[i][1] += t2tmp;
    torque[i][2] += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutDipoleLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCutDipoleLong::memory_usage();
  return bytes;
}
