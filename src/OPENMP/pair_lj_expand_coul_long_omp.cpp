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
   OpenMP-threaded variant of pair_lj_expand_coul_long.
   Follows the ThrOMP pattern established in pair_lj_cut_coul_long_omp.
   RESPA is disabled (respa_enable = 0) for simplicity.
------------------------------------------------------------------------- */

#include "pair_lj_expand_coul_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace EwaldConst;

namespace { using dbl3_t = struct { double x,y,z; }; }

/* ---------------------------------------------------------------------- */

PairLJExpandCoulLongOMP::PairLJExpandCoulLongOMP(LAMMPS *lmp) :
  PairLJExpandCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJExpandCoulLongOMP::compute(int eflag, int vflag)
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
void PairLJExpandCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj   = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int ** const firstneigh             = list->firstneigh;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    const double qtmp = q[i];
    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const int itype    = type[i];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    const double * _noalias const cutsqi    = cutsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const shifti    = shift[itype];
    const double * _noalias const lj1i      = lj1[itype];
    const double * _noalias const lj2i      = lj2[itype];
    const double * _noalias const lj3i      = lj3[itype];
    const double * _noalias const lj4i      = lj4[itype];
    const double * _noalias const offseti   = offset[itype];

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

        double forcecoul = 0.0;
        double prefactor = 0.0;
        double erfc_val  = 0.0;
        int itable = 0;
        double fraction = 0.0;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            const double r  = sqrt(rsq);
            const double grij = g_ewald * r;
            const double expm2 = exp(-grij*grij);
            const double t  = 1.0 / (1.0 + EWALD_P*grij);
            erfc_val = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc_val + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
            double table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        }

        double forcelj = 0.0;
        double r6inv   = 0.0;
        double rshift  = 0.0;
        if (rsq < cut_ljsqi[jtype]) {
          const double r = sqrt(rsq);
          rshift = r - shifti[jtype];
          const double rshiftsq   = rshift*rshift;
          const double rshift2inv = 1.0/rshiftsq;
          r6inv   = rshift2inv*rshift2inv*rshift2inv;
          forcelj = r6inv * (lj1i[jtype]*r6inv - lj2i[jtype]);
          forcelj *= factor_lj/rshift/r;
        }

        const double fpair = forcecoul*r2inv + forcelj;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        double evdwl = 0.0, ecoul = 0.0;
        if (EFLAG) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc_val;
            else {
              double table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          }
          if (rsq < cut_ljsqi[jtype]) {
            evdwl = r6inv*(lj3i[jtype]*r6inv-lj4i[jtype]) - offseti[jtype];
            evdwl *= factor_lj;
          }
        }

        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                                  evdwl, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJExpandCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJExpandCoulLong::memory_usage();
  return bytes;
}
