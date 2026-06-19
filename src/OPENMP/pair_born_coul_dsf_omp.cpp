// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_born_coul_dsf_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBornCoulDSFOMP::PairBornCoulDSFOMP(LAMMPS *lmp) :
  PairBornCoulDSF(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBornCoulDSFOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
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
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBornCoulDSFOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const double * _noalias const special_lj = force->special_lj;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;

  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double prefactor,erfcc,erfcd,arg;
  double rexp;
  double fxtmp,fytmp,fztmp;

  const int nlocal = atom->nlocal;
  int j,jj,jnum,jtype;

  evdwl = ecoul = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int    * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const ai = a[itype];
    const double * _noalias const rhoi = rho[itype];
    const double * _noalias const sigmai = sigma[itype];
    const double * _noalias const ci = c[itype];
    const double * _noalias const di = d[itype];
    const double * _noalias const rhoinvi = rhoinv[itype];
    const double * _noalias const born1i = born1[itype];
    const double * _noalias const born2i = born2[itype];
    const double * _noalias const born3i = born3[itype];
    const double * _noalias const offseti = offset[itype];

    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    // self coulombic energy

    if (EFLAG) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally_thr(this,i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0,thr);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j]/r;
          arg = alpha * r;
          erfcd = MathSpecial::expmsq(arg);
          erfcc = MathSpecial::my_erfcx(arg) * erfcd;
          forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +
                                   r*f_shift) * r;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        if (rsq < cut_ljsqi[jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          r = sqrt(rsq);
          rexp = exp((sigmai[jtype]-r)*rhoinvi[jtype]);
          forceborn = born1i[jtype]*r*rexp - born2i[jtype]*r6inv
            + born3i[jtype]*r2inv*r6inv;
        } else forceborn = 0.0;

        fpair = (forcecoul + factor_lj*forceborn) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsqi[jtype]) {
            evdwl = ai[jtype]*rexp - ci[jtype]*r6inv +
              di[jtype]*r6inv*r2inv - offseti[jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairBornCoulDSFOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBornCoulDSF::memory_usage();

  return bytes;
}
