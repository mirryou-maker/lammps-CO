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
   Contributing author: Hendrik Heenen (hendrik.heenen@mytum.de)
------------------------------------------------------------------------- */

#include "pair_buck_coul_long_cs.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

static constexpr double EWALD_F =  1.12837917;
static constexpr double EWALD_P =  9.95473818e-1;
static constexpr double B0      = -0.1335096380159268;
static constexpr double B1      = -2.57839507e-1;
static constexpr double B2      = -1.37203639e-1;
static constexpr double B3      = -8.88822059e-3;
static constexpr double B4      = -5.80844129e-3;
static constexpr double B5      =  1.14652755e-1;

static constexpr double EPSILON = 1.0e-20;
static constexpr double EPS_EWALD = 1.0e-6;
static constexpr double EPS_EWALD_SQR = 1.0e-12;

namespace {
struct dbl3_t { double x, y, z; };
}

/* ---------------------------------------------------------------------- */

PairBuckCoulLongCS::PairBuckCoulLongCS(LAMMPS *lmp) : PairBuckCoulLong(lmp)
{
  ewaldflag = pppmflag = 1;
  writedata = 1;
  single_enable = 0;
  ftable = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulLongCS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double rsq,r2inv,r6inv,forcecoul,forcebuck,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,u;
  double r,rexp;
  double fxtmp,fytmp,fztmp;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  const int inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int * const * const firstneigh = list->firstneigh;

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
    const double * _noalias const cutsqi    = cutsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const rhoinvi   = rhoinv[itype];
    const double * _noalias const buck1i    = buck1[itype];
    const double * _noalias const buck2i    = buck2[itype];
    fxtmp = fytmp = fztmp = 0.0;

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
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            prefactor = qqrd2e * qtmp*q[j];
            if (factor_coul < 1.0) {
              // When bonded parts are being calculated a minimal distance (EPS_EWALD)
              // has to be added to the prefactor and erfc in order to make the
              // used approximation functions for the Ewald correction valid
              grij = g_ewald * (r+EPS_EWALD);
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= (r+EPS_EWALD);
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - (1.0-factor_coul));
              // Additionally r2inv needs to be accordingly modified since the later
              // scaling of the overall force shall be consistent
              r2inv = 1.0/(rsq + EPS_EWALD_SQR);
            } else {
              grij = g_ewald * r;
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= r;
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsqi[jtype]) {
          r = sqrt(rsq);
          r6inv = r2inv*r2inv*r2inv;
          rexp = exp(-r*rhoinvi[jtype]);
          forcebuck = buck1i[jtype]*r*rexp - buck2i[jtype]*r6inv;
        } else forcebuck = 0.0;
        fpair = (forcecoul + factor_lj*forcebuck) * r2inv;

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
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsqi[jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
              offset[itype][jtype];
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
