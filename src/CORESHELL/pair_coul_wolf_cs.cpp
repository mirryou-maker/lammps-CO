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

#include "pair_coul_wolf_cs.h"

#include "atom.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

static constexpr double EPSILON = 1.0e-20;

/* ---------------------------------------------------------------------- */

PairCoulWolfCS::PairCoulWolfCS(LAMMPS *lmp) : PairCoulWolf( lmp )
{
   single_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulWolfCS::compute(int eflag, int vflag)
{
  int i, j, ii, jj, jnum;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double fxtmp, fytmp, fztmp;
  double rsq, forcecoul, factor_coul;
  double prefactor;
  double r;
  double erfcc, erfcd, v_sh, dvdrr, e_self, e_shift, f_shift, qisq;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const int newton_pair = force->newton_pair;
  const double qqrd2e = force->qqrd2e;

  // self and shifted coulombic energy

  e_self = v_sh = 0.0;
  e_shift = erfc(alf * cut_coul) / cut_coul;
  f_shift = -(e_shift + 2.0 * alf / MY_PIS * exp(-alf * alf * cut_coul * cut_coul)) / cut_coul;

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
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    qisq = qtmp * qtmp;
    e_self = -(e_shift / 2.0 + alf / MY_PIS) * qisq * qqrd2e;
    if (evflag) ev_tally(i, i, nlocal, 0, 0.0, e_self, 0.0, 0.0, 0.0, 0.0);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_coulsq) {
        rsq += EPSILON;
        // Add EPSILON for case: r = 0; Interaction must be removed
        // by special bond
        r = sqrt(rsq);
        prefactor = qqrd2e * qtmp * q[j] / r;
        erfcc = erfc(alf * r);
        erfcd = exp(-alf * alf * r * r);
        v_sh = (erfcc - e_shift * r) * prefactor;
        dvdrr = (erfcc / rsq + 2.0 * alf / MY_PIS * erfcd / r) + f_shift;
        forcecoul = dvdrr * rsq * prefactor;
        if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
        fpair = forcecoul / rsq;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (eflag) {
          ecoul = v_sh;
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        } else
          ecoul = 0.0;

        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

