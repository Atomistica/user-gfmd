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
#include "fix_contact_sphere.h"
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

const char syntax[] =
  "fix contact/sphere: Illegal fix command, syntax is 'fix ID group-ID "
  "contact/sphere <center-x> <center-y> <center-z> <radius> <epsilon> <sigma> "
  "<cutoff>'";

/* ---------------------------------------------------------------------- */

FixContactSphere::FixContactSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int iarg;
  char *endptr;
  char errstr[120];
  
  if (narg != 10)
    error->all(FLERR,syntax);

  cxvar = NULL;
  cyvar = NULL;
  czvar = NULL;

  if (!strcmp(arg[3], "CENTER")) {
    cx = domain->xprd_half;
  }
  else {
    if (strstr(arg[3], "v_") == arg[3]) {
      cxvar = strdup(arg[3]+2);
    }
    else {
      cx = strtod(arg[3], &endptr);
      if (endptr == arg[3])
        error->all(FLERR,syntax);
    }
  }

  if (!strcmp(arg[4], "CENTER")) {
    cy = domain->yprd_half;
  }
  else {
    if (strstr(arg[4], "v_") == arg[4]) {
      cyvar = strdup(arg[4]+2);
    }
    else {
      cy = strtod(arg[4], &endptr);
      if (endptr == arg[4])
        error->all(FLERR,syntax);
    }
  }

  if (!strcmp(arg[5], "CENTER")) {
    cz = domain->zprd_half;
  }
  else {
    if (strstr(arg[5], "v_") == arg[5]) {
      czvar = strdup(arg[5]+2);
    }
    else {
      cz = strtod(arg[5], &endptr);
      if (endptr == arg[5])
        error->all(FLERR,syntax);
    }
  }

  radius = strtod(arg[6], &endptr);
  if (endptr == arg[6])
    error->all(FLERR,syntax);
  epsilon = strtod(arg[7], &endptr);
  if (endptr == arg[7])
    error->all(FLERR,syntax);
  sigma = strtod(arg[8], &endptr);
  if (endptr == arg[8])
    error->all(FLERR,syntax);
  cutoff = strtod(arg[9], &endptr);
  if (endptr == arg[9])
    error->all(FLERR,syntax);

  // precomputed stuff

  radius_sq = radius*radius;
  radius_plus_cutoff_sq = (radius+cutoff)*(radius+cutoff);

  coeff1 = 48.0 * epsilon * pow(sigma,12.0);
  coeff2 = 24.0 * epsilon * pow(sigma,6.0);
  coeff3 = 4.0 * epsilon * pow(sigma,12.0);
  coeff4 = 4.0 * epsilon * pow(sigma,6.0);
  double r2inv = 1.0/(cutoff*cutoff);
  double r6inv = r2inv*r2inv*r2inv;
  offset = r6inv*(coeff3*r6inv - coeff4);

  // fix behavior

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // other stuff

  force_flag = 0;
  fsphere_loc[0] = fsphere_loc[1] = fsphere_loc[2] = fsphere_loc[3] = 
    fsphere_loc[4] = 0.0;
  fsphere[0] = fsphere[1] = fsphere[2] = fsphere[3] = fsphere[4] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixContactSphere::~FixContactSphere()
{
  if (cxvar)  free(cxvar);
  if (cyvar)  free(cyvar);
  if (czvar)  free(czvar);
}

/* ---------------------------------------------------------------------- */

int FixContactSphere::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  nevery = 1;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactSphere::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all(FLERR,"fix contact/sphere: RESPA not yet supported.");
}

/* ---------------------------------------------------------------------- */

void FixContactSphere::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactSphere::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactSphere::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (cxvar)  cx = input->variable->compute_equal(cxvar);
  if (cyvar)  cy = input->variable->compute_equal(cyvar);
  if (czvar)  cz = input->variable->compute_equal(czvar);

  fsphere_loc[0] = fsphere_loc[1] = fsphere_loc[2] = fsphere_loc[3] = 
    fsphere_loc[4] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double rx  = x[i][0] - cx;
      double ry  = x[i][1] - cy;
      double rz  = x[i][2] - cz;

      domain->minimum_image(rx, ry, rz);

      double rsq = rx*rx + ry*ry + rz*rz;

      if (rsq < radius_sq) {
        error->one(FLERR,"fix contact/sphere: atom inside sphere");
      }

      if (rsq < radius_plus_cutoff_sq) {
        double r = sqrt(rsq);
        double rinv = 1.0/(r - radius);
        double r2inv = rinv*rinv;
        double r6inv = r2inv*r2inv*r2inv;
        double df = r6inv*(coeff1*r6inv - coeff2) * rinv;
        double e = r6inv*(coeff3*r6inv - coeff4) - offset;
        double fx, fy, fz;

        fx = df*rx/r;
        fy = df*ry/r;
        fz = df*rz/r;
        
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        fsphere_loc[0] += e;
        fsphere_loc[1] -= fx;
        fsphere_loc[2] -= fy;
        fsphere_loc[3] -= fz;
        
        // count number of contacting atoms
        fsphere_loc[4] += 1.0;
      }
    }
  }

  force_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixContactSphere::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return the contact area
------------------------------------------------------------------------- */

double FixContactSphere::compute_scalar()
{
  if (force_flag == 0) {
    MPI_Allreduce(fsphere_loc, fsphere, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsphere[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixContactSphere::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(fsphere_loc, fsphere, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsphere[n+1];
}

