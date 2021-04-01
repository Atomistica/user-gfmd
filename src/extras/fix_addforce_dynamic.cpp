/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_addforce_dynamic.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddForceDynamic::FixAddForceDynamic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Not enough arguments for fix "
                            "addforce/dynamic command");

  if (!atom->gfmd_flag) {
    error->all(FLERR,"fix addforce/dynamic requires atom attributes gid and "
               "xeq; use gfmd atom style.");
  }

  energy_global_flag = 1;
  vector_flag = 1;
  size_vector = 12;
  global_freq = 1;
  extvector = 1;

  xvalue_ = atof(arg[3]);
  yvalue_ = atof(arg[4]);
  zvalue_ = atof(arg[5]);

  nx_ = atoi(arg[6]);
  ny_ = atoi(arg[7]);

  freq_ = atof(arg[8]);

  /* Reset all counters */

  nsteps_ = 0;

  vec_is_reduced_ = 0;
  for (int k = 0; k < 12; k++)  vec_[k] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixAddForceDynamic::~FixAddForceDynamic()
{
}

/* ---------------------------------------------------------------------- */

int FixAddForceDynamic::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::post_force(int vflag)
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double **x = atom->x;
  double **xeq = atom->xeq;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = update->dt*update->ntimestep;

  double xfc[3], xfs[3], vfc[3], vfs[3];
  for (int k = 0; k < 3; k++) {
    xfc[k] = 0.0; xfs[k] = 0.0; vfc[k] = 0.0; vfs[k] = 0.0;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double c = cos(2*M_PI*(xeq[i][0]*nx_/xprd+xeq[i][1]*ny_/yprd-freq_*t));
      double s = sin(2*M_PI*(xeq[i][0]*nx_/xprd+xeq[i][1]*ny_/yprd-freq_*t));

      for (int k = 0; k < 3; k++) {
        xfc[k] += c*x[i][k];
        xfs[k] += s*x[i][k];
        vfc[k] += c*v[i][k];
        vfs[k] += s*v[i][k];
      }

      f[i][0] += c*xvalue_;
      f[i][1] += c*yvalue_;
      f[i][2] += c*zvalue_;

    }
  }

  nsteps_++;

  vec_is_reduced_ = 0;
  for (int k = 0; k < 3; k++) {
    vec_[k] += xfc[k];
    vec_[k+3] += xfs[k];
    vec_[k+6] += vfc[k];
    vec_[k+9] += vfs[k];
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return the Fourier components of positions and velocities
------------------------------------------------------------------------- */

double FixAddForceDynamic::compute_vector(int n)
{
  // only sum across procs one time

  if (vec_is_reduced_ == 0) {
    MPI_Allreduce(vec_, reduced_vec_, 12, MPI_DOUBLE, MPI_SUM, world);
    if (nsteps_ > 0) {
      for (int k = 0; k < 12; k++)
        reduced_vec_[k] /= nsteps_;
    }
    vec_is_reduced_ = 1;
  }
  return reduced_vec_[n];
}
