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
   Contributing author: Ariel Lozano (arielzn@gmail.com)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include "pair_born_coul_dsf_cs.h"

#include "atom.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double EPSILON = 1.0e-20;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

/* ---------------------------------------------------------------------- */

PairBornCoulDSFCS::PairBornCoulDSFCS(LAMMPS *lmp) : PairBornCoulDSF(lmp)
{
  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBornCoulDSFCS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fxtmp,fytmp,fztmp;
  double r,rsq,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double prefactor,erfcc,erfcd,arg;
  double rexp;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const double * _noalias const special_coul = force->special_coul;
  const int newton_pair = force->newton_pair;
  const double qqrd2e = force->qqrd2e;

  const int inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
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
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

  // self coulombic energy
    if (eflag) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
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
        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j] / r;
          arg = alpha * r ;
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
        if (newton_pair || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (eflag) {
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

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

