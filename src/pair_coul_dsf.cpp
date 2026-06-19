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
   Contributing author: Trung Dac Nguyen (ORNL)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include "pair_coul_dsf.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using namespace MathConst;

namespace {
struct dbl3_t { double x, y, z; };
}

/* ---------------------------------------------------------------------- */

PairCoulDSF::PairCoulDSF(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairCoulDSF::~PairCoulDSF()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulDSF::compute(int eflag, int vflag)
{
  int i, j, ii, jj, jnum;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double r, rsq, forcecoul, factor_coul;
  double prefactor, erfcc, erfcd, t;
  double fxtmp, fytmp, fztmp;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
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
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];

    if (eflag) {
      double e_self = -(e_shift / 2.0 + alpha / MY_PIS) * qtmp * qtmp * qqrd2e;
      ev_tally(i, i, nlocal, 0, 0.0, e_self, 0.0, 0.0, 0.0, 0.0);
    }

    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_coulsq) {
        r = sqrt(rsq);
        prefactor = qqrd2e * qtmp * q[j] / r;
        erfcd = exp(-alpha * alpha * rsq);
        t = 1.0 / (1.0 + EWALD_P * alpha * r);
        erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

        forcecoul = prefactor * (erfcc / r + 2.0 * alpha / MY_PIS * erfcd + r * f_shift) * r;
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
          ecoul = prefactor * (erfcc - r * e_shift - rsq * f_shift);
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

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulDSF::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  alpha = utils::numeric(FLERR, arg[0], false, lmp);
  cut_coul = utils::numeric(FLERR, arg[1], false, lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulDSF::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulDSF::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style coul/dsf requires atom attribute q");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;
  double erfcc = erfc(alpha * cut_coul);
  double erfcd = exp(-alpha * alpha * cut_coul * cut_coul);
  f_shift = -(erfcc / cut_coulsq + 2.0 / MY_PIS * alpha * erfcd / cut_coul);
  e_shift = erfcc / cut_coul - f_shift * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulDSF::init_one(int /*i*/, int /*j*/)
{
  return cut_coul;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) { fwrite(&setflag[i][j], sizeof(int), 1, fp); }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDSF::write_restart_settings(FILE *fp)
{
  fwrite(&alpha, sizeof(double), 1, fp);
  fwrite(&cut_coul, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &alpha, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_coul, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_coul, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairCoulDSF::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                           double factor_coul, double /*factor_lj*/, double &fforce)
{
  double r, erfcc, erfcd, prefactor, t;
  double forcecoul, phicoul;

  forcecoul = phicoul = 0.0;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = force->qqrd2e * atom->q[i] * atom->q[j] / r;
    erfcd = exp(-alpha * alpha * rsq);
    t = 1.0 / (1.0 + EWALD_P * alpha * r);
    erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

    forcecoul = prefactor * (erfcc / r + 2.0 * alpha / MY_PIS * erfcd + r * f_shift) * r;
    if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

    phicoul = prefactor * (erfcc - r * e_shift - rsq * f_shift);
    if (factor_coul < 1.0) phicoul -= (1.0 - factor_coul) * prefactor;
  }

  fforce = forcecoul / rsq;

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulDSF::extract(const char *str, int &dim)
{
  if (strcmp(str, "cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return nullptr;
}
