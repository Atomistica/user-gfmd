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

/* ----------------------------------------------------------------------
   Contributing author: Lars Pastewka (JHU/Fh-IWM)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "respa.h"

#include "fix_fire.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixFIRE::FixFIRE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  minsteps = 10;
  incfac = 1.2;
  decfac = 0.5;
  mix_in = 0.1;
  mixdec = 0.99;
  max_dt = 0.001;

  if (narg == 10) {
    minsteps = force->inumeric(FLERR,arg[3]);
    incfac = force->numeric(FLERR,arg[4]);
    decfac = force->numeric(FLERR,arg[5]);
    mix_in = force->numeric(FLERR,arg[6]);
    mixdec = force->numeric(FLERR,arg[7]);
    max_dt = force->numeric(FLERR,arg[8]);
    xlimit = force->numeric(FLERR,arg[9]);
  }
  else {
    if (narg != 3) error->all(FLERR,"Illegal fix fire command. Parameters are: "
			      "fix ID group-ID fire [minsteps incfac decfac "
			      "min_in mixdec max_dt limit_dx]");
  }

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixFIRE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFIRE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  vlimitsq = (xlimit/dtv) * (xlimit/dtv); // TAS

  cut = minsteps;
  mix = mix_in;

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"respa not supported");
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixFIRE::initial_integrate(int vflag)
{
  double dtfm,vsq,scale; // vsq and scale needed to limit dx to help fire stability; see LAMMPS's min_fire

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]; // Limit dx
        if ((xlimit > 0) && (vsq > vlimitsq)) { 
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }


        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]; // Limit dx
        if ((xlimit > 0) && (vsq > vlimitsq)) {  
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        } 
 

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixFIRE::final_integrate()
{
  double dtfm,vsq,scale;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]; // Limit dx
        if ((xlimit > 0) && (vsq > vlimitsq)) {  
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        } 

      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]; // Limit dx
        if ((xlimit > 0) && (vsq > vlimitsq)) {
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }


      }
  }

  // Turn the global velocity a little bit  more along the global force...

  double vf1, vg_dot_vg1, Fg_dot_Fg1;
  vf1 = vg_dot_vg1 = Fg_dot_Fg1 = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vf1 += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
      vg_dot_vg1 += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      Fg_dot_Fg1 += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
    }
  }

  double vf, vg_dot_vg, Fg_dot_Fg;
  MPI_Allreduce(&vf1, &vf, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&vg_dot_vg1, &vg_dot_vg, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&Fg_dot_Fg1, &Fg_dot_Fg, 1, MPI_DOUBLE, MPI_SUM, world);

  //  printf("%i, vf1 = %e, vf = %e\n", comm->me, vf1, vf);
  //  printf("%i, vg_dot_vg1 = %f, vg = %f\n", comm->me, vg_dot_vg1, vg_dot_vg1);
  //  printf("%i, Fg = %f, Fg = %f\n", comm->me, Fg_dot_Fg1, Fg_dot_Fg1);

  double help = 0.0;
  if (Fg_dot_Fg > 0.0) {
    help = mix*sqrt(vg_dot_vg/Fg_dot_Fg);
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      v[i][0] = (1-mix)*v[i][0] + help*f[i][0];
      v[i][1] = (1-mix)*v[i][1] + help*f[i][1];
      v[i][2] = (1-mix)*v[i][2] + help*f[i][2];
    }
  }

  // Cut the velocities if the total power done by forces is negative

  if (vf < 0.0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	v[i][0] = 0.0;
	v[i][1] = 0.0;
	v[i][2] = 0.0;
      }
    }
    cut = minsteps;
    dtv = dtv*decfac;
    mix = mix_in;
  }
  else {
    // mixing is important only right after cut...
    if (cut < 0) {
      dtv = MIN(dtv*incfac, max_dt);
      mix = mix*mixdec;
    }
    else {
      cut = cut - 1;
    }
  }

  if (dtv != update->dt) {
    for (int i = 0; i < modify->nfix; i++) modify->fix[i]->reset_dt();
  }
}

/* ---------------------------------------------------------------------- */

void FixFIRE::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv);
}
