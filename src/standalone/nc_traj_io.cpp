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
   Contributing author: Lars Pastewka (Johns Hopkins University)
------------------------------------------------------------------------- */

#include <netcdf.h>

#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "nc_traj_io.h"

using namespace LAMMPS_NS;

enum{INT,DOUBLE};  // same as in dump_custom.cpp

const char NC_FRAME_STR[]         = "frame";
const char NC_SPATIAL_STR[]       = "spatial";
const char NC_ATOM_STR[]          = "atom";
const char NC_CELL_SPATIAL_STR[]  = "cell_spatial";
const char NC_CELL_ANGULAR_STR[]  = "cell_angular";
const char NC_LABEL_STR[]         = "label";

const char NC_TIME_STR[]          = "time";
const char NC_CELL_LENGTHS_STR[]  = "cell_lengths";
const char NC_CELL_ANGLES_STR[]   = "cell_angles";
const char NC_EPOT_STR[]          = "c_pe";
const char NC_EVDWL_STR[]         = "f_rigid";
const char NC_EEL0_STR[]          = "f_gfmd0";
const char NC_EEL_STR[]           = "f_gfmd";

const char NC_UNITS_STR[]         = "units";
const char NC_SCALE_FACTOR_STR[]  = "scale_factor";

const char NC_Z_STR[]             = "Z";
const char NC_COORDINATES_STR[]   = "coordinates";
const char NC_FORCES_STR[]        = "f_force";

const int MAX_DIMS           = 10;

/* ---------------------------------------------------------------------- */

#define NCERR(x) ncerr(x, __LINE__)

/* ---------------------------------------------------------------------- */

NCTrajIO::NCTrajIO(Pointers *lmp) : Pointers(lmp)
{
  double_precision = false;
}

/* ---------------------------------------------------------------------- */

NCTrajIO::~NCTrajIO()
{
  NCERR( nc_close(ncid) );
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::open_for_read(char *filename)
{
  // get total number of atoms
  ntotal = atom->nlocal;

  int dims[NC_MAX_VAR_DIMS];
  size_t index[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  double d[1];

  NCERR( nc_open(filename, NC_NOWRITE, &ncid) );
    
  // coordinates
  NCERR( nc_inq_varid(ncid, NC_COORDINATES_STR, &coordinates_var) );

  // forces
  NCERR( nc_inq_varid(ncid, NC_FORCES_STR, &forces_var) );
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::read_frame(int framei)
{
  size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  //  ptrdiff_t stride[NC_MAX_VAR_DIMS];

  start[0] = framei;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = ntotal;
  count[2] = 3;

  //  stride[0] = 1;
  //  stride[1] = 1;
  //  stride[2] = 3;

  NCERR( nc_get_vara_double(ncid, coordinates_var, start, count,
			    &atom->x[0][0]) );
  NCERR( nc_get_vara_double(ncid, forces_var, start, count,
			    &atom->f[0][0]) );
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::open_for_write(char *filename)
{
  // get total number of atoms
  ntotal = atom->nlocal;

  int dims[NC_MAX_VAR_DIMS];
  size_t index[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  double d[1];

  NCERR( nc_create(filename, NC_64BIT_OFFSET, &ncid) );
    
  // dimensions
  NCERR( nc_def_dim(ncid, NC_FRAME_STR, NC_UNLIMITED, &frame_dim) );
  NCERR( nc_def_dim(ncid, NC_SPATIAL_STR, 3, &spatial_dim) );
  NCERR( nc_def_dim(ncid, NC_ATOM_STR, ntotal, &atom_dim) );
  NCERR( nc_def_dim(ncid, NC_CELL_SPATIAL_STR, 3, &cell_spatial_dim) );
  NCERR( nc_def_dim(ncid, NC_CELL_ANGULAR_STR, 3, &cell_angular_dim) );
  NCERR( nc_def_dim(ncid, NC_LABEL_STR, 10, &label_dim) );

  // default variables
  dims[0] = spatial_dim;
  NCERR( nc_def_var(ncid, NC_SPATIAL_STR, NC_CHAR, 1, dims, &spatial_var) );
  NCERR( nc_def_var(ncid, NC_CELL_SPATIAL_STR, NC_CHAR, 1, dims,
		    &cell_spatial_var) );
  dims[0] = spatial_dim;
  dims[1] = label_dim;
  NCERR( nc_def_var(ncid, NC_CELL_ANGULAR_STR, NC_CHAR, 2, dims,
		    &cell_angular_var) );
  
  dims[0] = frame_dim;
  NCERR( nc_def_var(ncid, NC_TIME_STR, NC_DOUBLE, 1, dims, &time_var) );
  NCERR( nc_def_var(ncid, NC_EPOT_STR, NC_DOUBLE, 1, dims, &epot_var) );
  NCERR( nc_def_var(ncid, NC_EVDWL_STR, NC_DOUBLE, 1, dims, &evdwl_var) );
  NCERR( nc_def_var(ncid, NC_EEL0_STR, NC_DOUBLE, 1, dims, &eel0_var) );
  NCERR( nc_def_var(ncid, NC_EEL_STR, NC_DOUBLE, 1, dims, &eel_var) );
  dims[0] = frame_dim;
  dims[1] = cell_spatial_dim;
  NCERR( nc_def_var(ncid, NC_CELL_LENGTHS_STR, NC_DOUBLE, 2, dims,
		    &cell_lengths_var) );
  dims[0] = frame_dim;
  dims[1] = cell_angular_dim;
  NCERR( nc_def_var(ncid, NC_CELL_ANGLES_STR, NC_DOUBLE, 2, dims,
		    &cell_angles_var) );

  // variables specified in the input file
  dims[0] = frame_dim;
  dims[1] = atom_dim;
  dims[2] = spatial_dim;

  // Z
  NCERR( nc_def_var(ncid, NC_Z_STR, NC_INT, 1, dims+1, &Z_var) );

  // coordinates
  NCERR( nc_def_var(ncid, NC_COORDINATES_STR, NC_DOUBLE, 3, dims,
		    &coordinates_var) );

  // forces
  NCERR( nc_def_var(ncid, NC_FORCES_STR, NC_DOUBLE, 3, dims, &forces_var) );

  // attributes
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "Conventions",
			 5, "AMBER") );
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion",
			 3, "1.0") );
  
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "program",
			 7, "contact") );
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "programVersion",
			 3, "0.1") );

  // units
  NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
  			 2, "lj") );
  NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			 2, "lj") );
  
  NCERR( nc_put_att_text(ncid, cell_angles_var, NC_UNITS_STR,
			 6, "degree") );
  
  d[0] = 1.0;
  NCERR( nc_put_att_double(ncid, time_var, NC_SCALE_FACTOR_STR,
  			   NC_DOUBLE, 1, d) );
  d[0] = 1.0;
  NCERR( nc_put_att_double(ncid, cell_lengths_var, NC_SCALE_FACTOR_STR,
			   NC_DOUBLE, 1, d) );

  /*
   * Finished with definition
   */

  NCERR( nc_enddef(ncid) );

  /*
   * Write label variables
   */
  
  NCERR( nc_put_var_text(ncid, spatial_var, "xyz") );
  NCERR( nc_put_var_text(ncid, cell_spatial_var, "abc") );
  index[0] = 0;
  index[1] = 0;
  count[0] = 1;
  count[1] = 5;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "alpha") );
  index[0] = 1;
  count[1] = 4;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "beta") );
  index[0] = 2;
  count[1] = 5;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "gamma") );
    
  framei = 0;
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::write_header(double epot, double evdwl, double eel0, double eel)
{
  int *type = atom->type;

  size_t start[2];

  start[0] = framei;
  start[1] = 0;

  size_t count[2];
  double cell_lengths[3], cell_angles[3];

  double time = framei;
  if (domain->triclinic == 0) {
    cell_lengths[0] = domain->xprd;
    cell_lengths[1] = domain->yprd;
    cell_lengths[2] = domain->zprd;

    cell_angles[0] = 90;
    cell_angles[1] = 90;
    cell_angles[2] = 90;
  }
  else {
    error->all("NCTrajIO::write_header: Implement support for triclinic\n");
  }

  count[0] = 1;
  count[1] = 3;
  NCERR( nc_put_var1_double(ncid, time_var, start, &time) );
  NCERR( nc_put_vara_double(ncid, cell_lengths_var, start, count,
			    cell_lengths) );
  NCERR( nc_put_vara_double(ncid, cell_angles_var, start, count,
			    cell_angles) );
  NCERR( nc_put_var1_double(ncid, epot_var, start, &epot) );
  NCERR( nc_put_var1_double(ncid, evdwl_var, start, &evdwl) );
  NCERR( nc_put_var1_double(ncid, eel0_var, start, &eel0) );
  NCERR( nc_put_var1_double(ncid, eel_var, start, &eel) );

  // Z
  if (framei == 0) {
    if (!type) {
      count[0] = ntotal;
      type = new int[ntotal];
      for (int i = 0; i < ntotal; i++)
	type[i] = 2;
      NCERR( nc_put_vara_int(ncid, Z_var, start, count, type) );
      delete type;
    }
    else {
      NCERR( nc_put_vara_int(ncid, Z_var, start, count, type) );
    }
  }
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::write_data()
{
  size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  //  ptrdiff_t stride[NC_MAX_VAR_DIMS];

  start[0] = framei;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = ntotal;
  count[2] = 3;

  //  stride[0] = 1;
  //  stride[1] = 1;
  //  stride[2] = 3;

  NCERR( nc_put_vara_double(ncid, coordinates_var, start, count,
			    &atom->x[0][0]) );
  NCERR( nc_put_vara_double(ncid, forces_var, start, count,
			    &atom->f[0][0]) );

  NCERR( nc_sync(ncid) );
  framei++;
}

/* ---------------------------------------------------------------------- */

void NCTrajIO::ncerr(int err, int line)
{
  if (err != NC_NOERR) {
    char errstr[1024];
    sprintf(errstr, "NetCDF failed with error '%s' in line %i of %s.",
	    nc_strerror(err), line, __FILE__);
    error->one(errstr);
  }
}
