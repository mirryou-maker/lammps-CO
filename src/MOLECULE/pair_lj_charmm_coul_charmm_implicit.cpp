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

#include "pair_lj_charmm_coul_charmm_implicit.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

namespace {
using dbl3_t = struct { double x, y, z; };
}    // namespace

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmImplicit::PairLJCharmmCoulCharmmImplicit(LAMMPS *lmp) :
  PairLJCharmmCoulCharmm(lmp)
{
  implicit = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicit::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double fxtmp,fytmp,fztmp;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double philj,switch1,switch2;
  int * const *firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  const int inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;
    const double * _noalias const lj1i = lj1[itype];
    const double * _noalias const lj2i = lj2[itype];
    const double * _noalias const lj3i = lj3[itype];
    const double * _noalias const lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          forcecoul = 2.0 * qqrd2e * qtmp*q[j]*r2inv;
          if (rsq > cut_coul_innersq) {
            switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
              (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / denom_coul;
            switch2 = 12.0*rsq * (cut_coulsq-rsq) *
              (rsq-cut_coul_innersq) / denom_coul;
             forcecoul *= switch1 + 0.5*switch2;
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          jtype = type[j];
          forcelj = r6inv * (lj1i[jtype]*r6inv - lj2i[jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
            switch2 = 12.0*rsq * (cut_ljsq-rsq) *
              (rsq-cut_lj_innersq) / denom_lj;
            philj = r6inv * (lj3i[jtype]*r6inv - lj4i[jtype]);
            forcelj = forcelj*switch1 + philj*switch2;
          }
        } else forcelj = 0.0;

        fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

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
            ecoul = qqrd2e * qtmp*q[j]*r2inv;
            if (rsq > cut_coul_innersq) {
              switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
                (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) /
                denom_coul;
              ecoul *= switch1;
            }
            ecoul *= factor_coul;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq) {
            evdwl = r6inv*(lj3i[jtype]*r6inv-lj4i[jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
              evdwl *= switch1;
            }
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

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmImplicit::single(int i, int j,
                                              int itype, int jtype,
                                              double rsq,
                                              double factor_coul,
                                              double factor_lj,
                                              double &fforce)
{
  double r2inv,r6inv,switch1,switch2,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    forcecoul = 2.0 * force->qqrd2e * atom->q[i]*atom->q[j]*r2inv;
    if (rsq > cut_coul_innersq) {
      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
        (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / denom_coul;
      switch2 = 12.0*rsq * (cut_coulsq-rsq) *
        (rsq-cut_coul_innersq) / denom_coul;
       forcecoul *= switch1 + 0.5*switch2;
    }
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
        (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      switch2 = 12.0*rsq * (cut_ljsq-rsq) *
        (rsq-cut_lj_innersq) / denom_lj;
      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
      forcelj = forcelj*switch1 + philj*switch2;
    }
  } else forcelj = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*r2inv;
    if (rsq > cut_coul_innersq) {
      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
        (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) /
        denom_coul;
      phicoul *= switch1;
    }
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
        (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      philj *= switch1;
    }
    eng += factor_lj*philj;
  }

  return eng;
}
