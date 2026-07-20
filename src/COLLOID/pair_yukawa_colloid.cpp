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
   Contributing authors: Randy Schunk (Sandia)
------------------------------------------------------------------------- */

#include "pair_yukawa_colloid.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

/* ---------------------------------------------------------------------- */

PairYukawaColloid::PairYukawaColloid(LAMMPS *lmp) : PairYukawa(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

void PairYukawaColloid::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,radi,radj;
  double fxtmp,fytmp,fztmp;
  double rsq,r,rinv,screening,forceyukawa,factor;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const radius = atom->radius;
  const int * _noalias const type = atom->type;
  int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    radi = radius[i];
    const double * _noalias const cutsqi   = cutsq[itype];
    const double * _noalias const ai       = a[itype];
    const double * _noalias const offseti  = offset[itype];
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsqi[jtype]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        screening = exp(-kappa*(r-(radi+radj)));
        forceyukawa = ai[jtype] * screening;

        fpair = factor*forceyukawa * rinv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (eflag) {
          evdwl = ai[jtype]/kappa * screening - offseti[jtype];
          evdwl *= factor;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairYukawaColloid::init_style()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Pair yukawa/colloid requires atom attribute radius");

  neighbor->add_request(this);

  // require that atom radii are identical within each type

  for (int i = 1; i <= atom->ntypes; i++)
    if (!atom->radius_consistency(i,rad[i]))
      error->all(FLERR,"Pair yukawa/colloid requires atoms with same type have same radius");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairYukawaColloid::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a[i][j] = mix_energy(a[i][i],a[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag && (kappa != 0.0)) {
    double screening = exp(-kappa * (cut[i][j] - (rad[i]+rad[j])));
    offset[i][j] = a[i][j]/kappa * screening;
  } else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairYukawaColloid::single(int /*i*/, int /*j*/, int itype, int jtype,
                                 double rsq,
                                 double /*factor_coul*/, double factor_lj,
                                 double &fforce)
{
  double r,rinv,screening,forceyukawa,phi;

  r = sqrt(rsq);
  rinv = 1.0/r;
  screening = exp(-kappa*(r-(rad[itype]+rad[jtype])));
  forceyukawa = a[itype][jtype] * screening;
  fforce = factor_lj*forceyukawa * rinv;

  phi = a[itype][jtype]/kappa * screening  - offset[itype][jtype];
  return factor_lj*phi;
}
