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

#include "pair_morse_soft_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMorseSoftOMP::PairMorseSoftOMP(LAMMPS *lmp) :
  PairMorseSoft(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairMorseSoftOMP::compute(int eflag, int vflag)
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

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairMorseSoftOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const double * _noalias const special_lj = force->special_lj;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  double xtmp,ytmp,ztmp,delx,dely,delz,fxtmp,fytmp,fztmp;
  double rsq,r,dr,dexp,dexp2,dexp3,factor_lj,evdwl,fpair;
  double ea,phi,V0,iea2;
  double D,a,x0,l,B,s1,llf;

  const int nlocal = atom->nlocal;
  int j,jj,jnum,jtype;

  evdwl = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int    * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const d0i = d0[itype];
    const double * _noalias const alphai = alpha[itype];
    const double * _noalias const r0i = r0[itype];
    const double * _noalias const lambdai = lambda[itype];

    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r = sqrt(rsq);
        dr = r - r0i[jtype];

        D = d0i[jtype];
        a = alphai[jtype];
        x0 = r0i[jtype];
        dexp = exp(-a * dr);
        dexp2 = dexp * dexp;
        dexp3 = dexp2 * dexp;

        l = lambdai[jtype];

        ea = exp(a * x0);
        iea2 = exp(-2. * a * x0);

        V0 = D * dexp * (dexp - 2.0);
        B = -2.0 * D * iea2 * (ea - 1.0) / 3.0;

        if (l >= shift_range) {
          s1 = (l - 1.0) / (shift_range - 1.0);
          phi = V0 + B * dexp3 * s1;

          // Force computation:
          fpair = 3.0 * a * B * dexp3 * s1 + 2.0 * a * D * (dexp2 - dexp);
          fpair /= r;
        } else {
          llf = MathSpecial::powint(l / shift_range, nlambda);
          phi = V0 + B * dexp3;
          phi *= llf;

          // Force computation:
          if (r == 0.0) {
            fpair = 0.0;
          } else {
            fpair = 3.0 * a * B * dexp3 + 2.0 * a * D * (dexp2 - dexp);
            fpair *= llf / r;
          }
        }

        fpair *= factor_lj;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          evdwl = phi * factor_lj;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairMorseSoftOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMorseSoft::memory_usage();

  return bytes;
}
