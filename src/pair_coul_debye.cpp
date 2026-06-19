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

#include "pair_coul_debye.h"

#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "error.h"


using namespace LAMMPS_NS;

namespace {
struct dbl3_t { double x, y, z; };
}

/* ---------------------------------------------------------------------- */

PairCoulDebye::PairCoulDebye(LAMMPS *lmp) : PairCoulCut(lmp)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

void PairCoulDebye::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,r2inv,r,rinv,forcecoul,factor_coul,screening;
  double fxtmp,fytmp,fztmp;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
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
    itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const scalei = scale[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        rinv = 1.0/r;
        screening = exp(-kappa*r);
        forcecoul = qqrd2e * scalei[jtype] *
          qtmp*q[j] * screening * (kappa + rinv);
        fpair = factor_coul*forcecoul * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (eflag) ecoul = factor_coul * qqrd2e *
            scalei[jtype] * qtmp*q[j] * rinv * screening;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulDebye::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  kappa = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDebye::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDebye::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&kappa,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulDebye::single(int i, int j, int /*itype*/, int /*jtype*/,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double r2inv,r,rinv,forcecoul,phicoul,screening;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rinv = 1.0/r;
  screening = exp(-kappa*r);
  forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
    screening * (kappa + rinv);
  fforce = factor_coul*forcecoul * r2inv;

  phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening;
  return factor_coul*phicoul;
}

/* ---------------------------------------------------------------------- */

void PairCoulDebye::born_matrix(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                            double factor_coul, double /*factor_lj*/, double &dupair,
                            double &du2pair)
{
  double r, rinv, r2inv, r3inv, screening;
  double du_coul, du2_coul;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  r = sqrt(rsq);
  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  r3inv = r2inv * rinv;
  screening = exp(-kappa*r);

  // Reminder: qqrd2e converts  q^2/r to energy w/ dielectric constant
  du_coul = -qqrd2e * q[i] * q[j] * r2inv * (1 + kappa * r) *  screening;
  du2_coul = qqrd2e * q[i] * q[j] * r3inv * (2 + 2 * kappa * r + kappa * kappa * rsq) * screening;

  dupair = factor_coul * du_coul;
  du2pair = factor_coul * du2_coul;
}
