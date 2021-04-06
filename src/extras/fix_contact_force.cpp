/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

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

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <netcdfcpp.h>

#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include "fix_contact_force.h"

using namespace LAMMPS_NS;

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

const char syntax[] =
  "fix contact/force: Illegal fix command, syntax is 'fix ID group-ID "
  "contact/force <file>'";

#define MODULO(i, n)  ( (i) >= 0 ? (i) % (n) : (n+i) % (n) )

/* ---------------------------------------------------------------------- */

FixContactForce::FixContactForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4)
    error->all(syntax);

  int iarg = 3;

  // read forces

  read_forces(arg[iarg++]);

  // fix behavior

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // other stuff

  force_flag = 0;
  force_loc[0] = force_loc[1] = force_loc[2] = force_loc[3] = 0.0;
  force[0] = force[1] = force[2] = force[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixContactForce::~FixContactForce()
{
  memory->destroy_2d_double_array(forces);
}

/* ---------------------------------------------------------------------- */

int FixContactForce::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= THERMO_ENERGY;
  mask |= END_OF_STEP;
  nevery = 1;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactForce::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all("fix contact/map: RESPA not yet supported.");
}

/* ---------------------------------------------------------------------- */

void FixContactForce::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactForce::min_setup(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactForce::pre_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  force_loc[0] = force_loc[1] = force_loc[2] = force_loc[3] = 0.0;
  force_flag = 0;

  for (int i = 0; i < nlocal; i++) {
    int j = tag[i];

    f[i][0] += forces[j][0];
    f[i][1] += forces[j][1];
    f[i][2] += forces[j][2];

    force_loc[0] += forces[j][0];
    force_loc[1] += forces[j][1];
    force_loc[2] += forces[j][2];

    force_loc[3] += forces[j][0]*x[i][0] + forces[j][1]*x[i][1] +
      forces[j][2]*x[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixContactForce::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ----------------------------------------------------------------------
   return the contact area
------------------------------------------------------------------------- */

double FixContactForce::compute_scalar()
{
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(force_loc, force, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return force[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixContactForce::compute_vector(int n)
{
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(force_loc, force, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return force[n+1];
}

/* ----------------------------------------------------------------------
   parse file name into file and frame number
------------------------------------------------------------------------- */

int parse_fn(char *fn)
{
  int frame = 0;

  char *at = strchr(fn, '@');
  if (at) {
    frame = atoi(at+1);
    *at = 0;
  }

  return frame;
}

/* ----------------------------------------------------------------------
   read force from NetCDF file
------------------------------------------------------------------------- */

void FixContactForce::read_forces(char *in_fn)
{
  char *fn = strdup(in_fn);
  int frame;
  frame = parse_fn(fn);

  NcFile *nc = new NcFile(fn);

  int nlocal = atom->nlocal, nall;
  MPI_Allreduce(&nall, &nlocal, 1, MPI_INT, MPI_SUM, world);

  int nc_nframes = nc->get_dim("frame")->size();
  int nc_nall = nc->get_dim("atom")->size();

  if (frame >= nc_nframes) {
    char errstr[1024];
    sprintf(errstr, "fix contact/force: Given frame (%i) larger than total "
	    "number of frames in .nc file (%i).", frame, nc_nframes);
  }

  if (nc_nall != nall) {
    char errstr[1024];
    sprintf(errstr, "fix contact/force: Number of atoms (%i) differs from "
	    "number of atoms in .nc file (%i).", nall, nc_nall);
    error->all(errstr);
  }

  NcVar *f_var = nc->get_var("f_forces");

  forces = memory->create_2d_double_array(nall, 3, "FicContactForce::forces");
  if (!f_var->set_cur(frame))
    error->one("fix contact/force: set_cur failed.");
  if (!f_var->get(forces[0], 1, nall, 3))
    error->all("fix contact/force: get failed.");

  delete f_var;
  delete nc;

  free(fn);
}

