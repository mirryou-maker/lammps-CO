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
   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu
------------------------------------------------------------------------- */

#include "pair_mdpd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace {
  struct dbl3_t { double x,y,z; };
}

static constexpr double EPSILON = 1.0e-10;

static const char cite_pair_mdpd[] =
  "pair mdpd command: https://doi.org/10.1063/1.4812366\n\n"
  "@Article{ZLi2013_POF,\n"
  " author = {Li, Z. and Hu, G. H. and Wang, Z. L. and Ma Y. B. and Zhou, Z. W.},\n"
  " title = {Three Dimensional Flow Structures in a Moving Droplet on Substrate: a Dissipative Particle Dynamics Study},\n"
  " journal = {Physics of Fluids},\n"
  " year = {2013},\n"
  " volume = {25},\n"
  " number = {7},\n"
  " pages = {072103}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairMDPD::PairMDPD(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_mdpd);

  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairMDPD::~PairMDPD()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_r);
    memory->destroy(A_att);
    memory->destroy(B_rep);
    memory->destroy(gamma);
    memory->destroy(sigma);
  }
  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairMDPD::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wc,wc_r, wr,randnum,factor_dpd;
  int **firstneigh;
  double rhoi, rhoj;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const auto * _noalias const vel = (dbl3_t *) atom->v[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  double *rho= atom->rho;
  const int * _noalias const type = atom->type;
  int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    vxtmp = vel[i].x;
    vytmp = vel[i].y;
    vztmp = vel[i].z;
    itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    rhoi = rho[i];
    const double * _noalias const cutsqi  = cutsq[itype];
    const double * _noalias const cuti    = cut[itype];
    const double * _noalias const cut_ri  = cut_r[itype];
    const double * _noalias const A_atti  = A_att[itype];
    const double * _noalias const B_repi  = B_rep[itype];
    const double * _noalias const gammai  = gamma[itype];
    const double * _noalias const sigmai  = sigma[itype];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in MDPD systems
        rinv = 1.0/r;
        delvx = vxtmp - vel[j].x;
        delvy = vytmp - vel[j].y;
        delvz = vztmp - vel[j].z;
        dot = delx*delvx + dely*delvy + delz*delvz;

        wc = 1.0 - r/cuti[jtype];
        wc_r = 1.0 - r/cut_ri[jtype];
        wc_r = MAX(wc_r,0.0);
        wr = wc;

        rhoj = rho[j];
        randnum = random->gaussian();

        // conservative force = A_att * wc + B_rep*(rhoi+rhoj)*wc_r
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr * rnd * dtinvsqrt;

        fpair = A_atti[jtype]*wc + B_repi[jtype]*(rhoi+rhoj)*wc_r;
        fpair -= gammai[jtype]*wr*wr*dot*rinv;
        fpair += sigmai[jtype]*wr*randnum*dtinvsqrt;
        fpair *= factor_dpd*rinv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*A_atti[jtype]*cuti[jtype] * wr*wr + 0.5*B_repi[jtype]*cut_ri[jtype]*(rhoi+rhoj)*wc_r*wc_r;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
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

void PairMDPD::allocate()
{
  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_r,n+1,n+1,"pair:cut_r");
  memory->create(A_att,n+1,n+1,"pair:A_att");
  memory->create(B_rep,n+1,n+1,"pair:B_rep");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMDPD::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed
  // create a positive seed based on the system clock, if requested.

  if (seed <= 0) {
    constexpr double LARGE_NUM = 2<<30;
    seed = int(fmod(platform::walltime() * LARGE_NUM, LARGE_NUM)) + 1;
  }

  delete random;
  random = new RanMars(lmp,(seed + comm->me) % 900000000);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMDPD::coeff(int narg, char **arg)
{
  if (narg != 7 ) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double A_one = utils::numeric(FLERR,arg[2],false,lmp);
  double B_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[4],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[5],false,lmp);
  double cut_two = utils::numeric(FLERR,arg[6],false,lmp);

  if (cut_one < cut_two) error->all(FLERR, "Value for cutA should be larger than cutB.");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      A_att[i][j] = A_one;
      B_rep[i][j] = B_one;
      gamma[i][j] = gamma_one;
      cut[i][j] = cut_one;
      cut_r[i][j] = cut_two;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMDPD::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair mdpd requires ghost atoms store velocity");

  if (!atom->rho_flag)
    error->all(FLERR,"Pair style mdpd requires atom attribute rho");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR, "Pair mdpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status:\n" + Info::get_pair_coeff_status(lmp));

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut[j][i] = cut[i][j];
  cut_r[j][i] = cut_r[i][j];
  A_att[j][i] = A_att[i][j];
  B_rep[j][i] = B_rep[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&A_att[i][j],sizeof(double),1,fp);
        fwrite(&B_rep[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&cut_r[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&A_att[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&B_rep[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_r[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&A_att[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&B_rep[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_r[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&temperature,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMDPD::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,A_att[i][i],B_rep[i][i],gamma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMDPD::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,A_att[i][j],B_rep[i][j],gamma[i][j],cut[i][j],cut_r[i][j]);
}

