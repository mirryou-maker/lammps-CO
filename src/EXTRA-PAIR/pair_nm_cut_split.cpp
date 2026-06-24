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
   Contributing Author: Julien Devemy (ICCF), Robert S. Hoy (USF), Joseph D. Dietz (USF)
------------------------------------------------------------------------- */

#include "pair_nm_cut_split.h"

#include "atom.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathSpecial::powint;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

/* ---------------------------------------------------------------------- */
PairNMCutSplit::PairNMCutSplit(LAMMPS *lmp) : PairNMCut(lmp)
{
  writedata = 1;
}

void PairNMCutSplit::compute(int eflag, int vflag)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double fxtmp, fytmp, fztmp;
  double rsq, r2inv, factor_lj;
  double r, forcenm, rminv, rninv;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const int newton_pair = force->newton_pair;

  const int inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const r0i = r0[itype];
    const double * _noalias const e0i = e0[itype];
    const double * _noalias const e0nmi = e0nm[itype];
    const double * _noalias const nmi = nm[itype];
    const double * _noalias const r0ni = r0n[itype];
    const double * _noalias const r0mi = r0m[itype];
    const double * _noalias const nni = nn[itype];
    const double * _noalias const mmi = mm[itype];
    const double * _noalias const offseti = offset[itype];

    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r2inv = 1.0 / rsq;
        r = sqrt(rsq);

        // r < r0 --> use generalized LJ
        if (rsq < r0i[jtype] * r0i[jtype]) {
          rminv = pow(r2inv, mmi[jtype] / 2.0);
          rninv = pow(r2inv, nni[jtype] / 2.0);
          // r^n = 1/rninv, r^m = 1/rminv (eliminates 2 pow calls)
          forcenm = e0nmi[jtype] * nmi[jtype] *
            (r0ni[jtype] * rninv - r0mi[jtype] * rminv);
        }
        // r > r0 --> use standard LJ (m = 6 n = 12)
        else forcenm = (e0i[jtype] / 6.0) * 72.0 * (4.0 / powint(r, 12) - 2.0 / powint(r, 6));

        fpair = factor_lj * forcenm * r2inv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (eflag) {
          // r < r0 --> use generalized LJ
          if (rsq < r0i[jtype] * r0i[jtype]) {
            // rminv/rninv already computed in force block above
            evdwl = e0nmi[jtype] * (mmi[jtype] * r0ni[jtype] * rninv -
                                    nni[jtype] * r0mi[jtype] * rminv) -
              offseti[jtype];
          }
          // r > r0 --> use standard LJ (m = 6 n = 12)
          else evdwl = (e0i[jtype] / 6.0) * (24.0 * powint(r2inv, 6) - 24.0 * pow(r2inv, 3));
          evdwl *= factor_lj;
        }
        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairNMCutSplit::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq, double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv,rminv,rninv,r,forcenm,phinm;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rminv = pow(r2inv,mm[itype][jtype]/2.0);
  rninv = pow(r2inv,nn[itype][jtype]/2.0);
  // r < r_0, use generalized LJ
  if (rsq < r0[itype][jtype]*r0[itype][jtype]) {  // note the addition of the r0 factor
     forcenm = e0nm[itype][jtype]*nm[itype][jtype]*
      (r0n[itype][jtype]*rninv - r0m[itype][jtype]*rminv);
      phinm = e0nm[itype][jtype]*(mm[itype][jtype]*r0n[itype][jtype]*rninv
      -nn[itype][jtype]*r0m[itype][jtype]*rminv)-offset[itype][jtype];

  }
  // r > r_0 --> use standard LJ (m = 6 n = 12)
  else {
    forcenm = (e0[itype][jtype]/6.0)*72.0*(4.0/powint(r,12)-2.0/powint(r,6));
    phinm = (e0[itype][jtype]/6.0)*(24.0*powint(r2inv,6)-24.0*powint(r2inv,3));
  }

  fforce = factor_lj*forcenm*r2inv;
  return factor_lj*phinm;
}

