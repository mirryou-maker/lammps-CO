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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "pair_gran_hooke.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "fix.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

/* ---------------------------------------------------------------------- */

PairGranHooke::PairGranHooke(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
  no_virial_fdotr_compute = 0;
  history = 0;
}

/* ---------------------------------------------------------------------- */

void PairGranHooke::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,factor_lj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,ft,fs1,fs2,fs3;
  ev_init(eflag,vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const auto * _noalias const v = (dbl3_t *) atom->v[0];
  auto * _noalias const f = (dbl3_t *) atom->f[0];
  const auto * _noalias const omega = (dbl3_t *) atom->omega[0];
  auto * _noalias const torque = (dbl3_t *) atom->torque[0];
  const double * _noalias const radius = atom->radius;
  const double * _noalias const rmass = atom->rmass;
  const int * _noalias const mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double * _noalias const special_lj = force->special_lj;

  int inum = list->inum;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    radi = radius[i];
    const int * _noalias const jlist = firstneigh[i];
    int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;
    double torxtmp = 0.0, torytmp = 0.0, torztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq < radsum*radsum) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i].x - v[j].x;
        vr2 = v[i].y - v[j].y;
        vr3 = v[i].z - v[j].z;

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi*omega[i].x + radj*omega[j].x) * rinv;
        wr2 = (radi*omega[i].y + radj*omega[j].y) * rinv;
        wr3 = (radi*omega[i].z + radj*omega[j].z) * rinv;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi*mj / (mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        // normal forces = Hookian contact + normal velocity damping

        damp = meff*gamman*vnnr*rsqinv;
        ccel = kn*(radsum-r)*rinv - damp;
        if (limit_damping && (ccel < 0.0)) ccel = 0.0;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // force normalization

        fn = xmu * fabs(ccel*r);
        fs = meff*gammat*vrel;
        if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
        else ft = 0.0;

        // tangential force due to tangential velocity damping

        fs1 = -ft*vtr1;
        fs2 = -ft*vtr2;
        fs3 = -ft*vtr3;

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        fx *= factor_lj;
        fy *= factor_lj;
        fz *= factor_lj;
        fxtmp += fx;
        fytmp += fy;
        fztmp += fz;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);
        tor1 *= factor_lj;
        tor2 *= factor_lj;
        tor3 *= factor_lj;
        torxtmp -= radi*tor1;
        torytmp -= radi*tor2;
        torztmp -= radi*tor3;

        if (newton_pair || j < nlocal) {
          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
          torque[j].x -= radj*tor1;
          torque[j].y -= radj*tor2;
          torque[j].z -= radj*tor3;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    torque[i].x += torxtmp;
    torque[i].y += torytmp;
    torque[i].z += torztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairGranHooke::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                             double /*factor_coul*/, double /*factor_lj*/,
                             double &fforce)
{
  double radi,radj,radsum,r,rinv,rsqinv;
  double delx,dely,delz;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel;
  double fn,fs,ft;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  // zero out forces if caller requests non-touching pair outside cutoff

  if (rsq >= radsum*radsum) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  const auto * const v_arr = (dbl3_t *) atom->v[0];
  vr1 = v_arr[i].x - v_arr[j].x;
  vr2 = v_arr[i].y - v_arr[j].y;
  vr3 = v_arr[i].z - v_arr[j].z;

  // normal component

  const auto * const x_arr = (dbl3_t *) atom->x[0];
  delx = x_arr[i].x - x_arr[j].x;
  dely = x_arr[i].y - x_arr[j].y;
  delz = x_arr[i].z - x_arr[j].z;

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  const auto * const omega_arr = (dbl3_t *) atom->omega[0];
  wr1 = (radi*omega_arr[i].x + radj*omega_arr[j].x) * rinv;
  wr2 = (radi*omega_arr[i].y + radj*omega_arr[j].y) * rinv;
  wr3 = (radi*omega_arr[i].z + radj*omega_arr[j].z) * rinv;

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: ensure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radsum-r)*rinv - damp;
  if (limit_damping && (ccel < 0.0)) ccel = 0.0;

  // relative velocities

  vtr1 = vt1 - (delz*wr2-dely*wr3);
  vtr2 = vt2 - (delx*wr3-delz*wr1);
  vtr3 = vt3 - (dely*wr1-delx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = meff*gammat*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // set force and return no energy

  fforce = ccel;


  // set single_extra quantities

  svector[0] = -ft*vtr1;
  svector[1] = -ft*vtr2;
  svector[2] = -ft*vtr3;
  svector[3] = sqrt(svector[0]*svector[0] +
                    svector[1]*svector[1] +
                    svector[2]*svector[2]);
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;

  return 0.0;
}
