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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "displace_xeq.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "irregular.h"
#include "group.h"
#include "random_park.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{MOVE_BY,MOVE_TO,RAMP};

/* ---------------------------------------------------------------------- */

DisplaceXEq::DisplaceXEq(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void DisplaceXEq::command(int narg, char **arg)
{
  int i;

  if (!atom->gfmd_flag) {
    error->all(FLERR,"displace_xeq command requires atom attributes gid and xeq; "
               "use gfmd atom style.");
  }

  if (domain->box_exist == 0) 
    error->all(FLERR,"displace_xeq command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace_xeq command");
  if (modify->nfix_restart_peratom) 
    error->all(FLERR,"Cannot displace_xeq after "
	       "reading restart file with per-atom info");

  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");

  // group and style

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_xeq group ID");
  int groupbit = group->bitmask[igroup];

  int style;
  if (strcmp(arg[1],"move_by") == 0) style = MOVE_BY;
  else if (strcmp(arg[1],"move_to") == 0) style = MOVE_TO;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else error->all(FLERR,"Use either move_by, move_to or ramp.");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE_BY || style == MOVE_TO) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all(FLERR,"Use of displace_xeq with undefined lattice.");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector

  if (style == MOVE_BY) {

    double delx = xscale*atof(arg[2]);
    double dely = yscale*atof(arg[3]);
    double delz = zscale*atof(arg[4]);

    double **xeq = atom->xeq;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;

    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;

    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
	xeq[i][0] += delx;
	xeq[i][1] += dely;
	xeq[i][2] += delz;
	sumx += xeq[i][0];
	sumy += xeq[i][1];
	sumz += xeq[i][2];
      }
    }

    if (comm->me == 0 && screen) fprintf(screen,"New avg xeq is %f %f %f\n",
					 sumx/nall, sumy/nall, sumz/nall);
  }

  // move atoms to z-coordinate

  if (style == MOVE_TO) {

    double delz = zscale*atof(arg[4]);

    double **xeq = atom->xeq;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;

    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
	xeq[i][2] = delz;
      }
    }
  }

  // move atoms in ramped fashion
    
  if (style == RAMP) {

    int d_dim;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all(FLERR,"Illegal displace_xeq ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*atof(arg[3]);
      d_hi = xscale*atof(arg[4]);
    } else if (d_dim == 1) {
      d_lo = yscale*atof(arg[3]);
      d_hi = yscale*atof(arg[4]);
    } else if (d_dim == 2) {
      d_lo = zscale*atof(arg[3]);
      d_hi = zscale*atof(arg[4]);
    }

    int coord_dim;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all(FLERR,"Illegal displace_xeq ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*atof(arg[6]);
      coord_hi = xscale*atof(arg[7]);
    } else if (coord_dim == 1) {
      coord_lo = yscale*atof(arg[6]);
      coord_hi = yscale*atof(arg[7]);
    } else if (coord_dim == 2) {
      coord_lo = zscale*atof(arg[6]);
      coord_hi = zscale*atof(arg[7]);
    }

    double **xeq = atom->xeq;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;

    double fraction,dramp;

    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
	fraction = (xeq[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
	fraction = MAX(fraction,0.0);
	fraction = MIN(fraction,1.0);
	dramp = d_lo + fraction*(d_hi - d_lo);
	xeq[i][d_dim] += dramp;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_xeq input line 
------------------------------------------------------------------------- */

void DisplaceXEq::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace_xeq command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace_xeq command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace_xeq command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace_xeq command");
  }
}
