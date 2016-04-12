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

// TAS Modifying fix_addforce/dynamic and fix_gravity and fix_addforce
FixAddForceDynamic::FixAddForceDynamic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Not enough arguments fix id addforce/dynamic "
			    "acontact shiftx shifty (all nn units) nndist r_overEstar (<-not used in uniform case)");

  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  size_vector = 3;
  extscalar = 1;
  extvector = 1;

  acontact_nn_ = atof(arg[3]);
  shiftx_nn_ = atof(arg[4]);
  shifty_nn_ = atof(arg[5]);
  nndist_ = atof(arg[6]);
  r_overEstar = atof(arg[7]);

  /* Reset all counters */

  //eflag = 0;
  //egrav = 0.0;

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

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
  mask |= THERMO_ENERGY;
  mask |= MIN_POST_FORCE; // Not in gravity.cpp but I want it // TAS
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::min_setup(int vflag)
{
  //if (logfile) fprintf(logfile, "ADF calling post force from min_setup\n");
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::setup(int vflag) // TAS added 4.17.14 to see if can dump fix forces
{
  //if (logfile) fprintf(logfile, "ADF calling post force from setup\n");
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    if (logfile) fprintf(logfile, " adf not defined for respa.  exit. \n");
    exit(999);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::post_force(int vflag)
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  //printf("** addforce dynamic is adding force \n"); // TAS
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  double cx, cy, sx, sy;
  double envelope, rad_sq;

  // Define envelope (contacting region)
  // center of box
  double xcenter_sigma = xprd / 2;
  double ycenter_sigma = yprd / 2;  
  double acontact_sigma = nndist_ * acontact_nn_;
  
  // Define contact potential
  double radspersigma_x = 2.0*M_PI / nndist_;
  double radspersigma_y = 2.0*M_PI / nndist_;
  double border_dist = 0.5 * nndist_;

  if ((acontact_sigma+border_dist) > xprd/2) {
    if (logfile) fprintf(logfile, 
      "contact diam %f larger than system %f. exit.\n", 2*acontact_sigma, xprd);
    exit(999); // and this even ignores possible gfmd grid deformation
  }
  if ((acontact_sigma+border_dist) > xprd/4) {
    if (logfile) fprintf(logfile, 
      "Warning contact diam larger > system/2 so nonperiodic inappropriate \n");
  }

  // Put many calculations in variables to avoid calculating each loop
  double aradsq = pow(acontact_sigma, 2.0);
  double arad_border_sq = pow(acontact_sigma + border_dist, 2.0);

  // Only for othogonal boxes

  //printf("In post_force freq... %f (wavelengths/sys), rads per sigma %f\n", spatialfreqx, radspersigma_x);
  //printf(" huh? radspersigma_x %f radspersigma_y %f xshift_rads %f yshift_rads %f \n", radspersigma_x, radspersigma_y, xshift_rads, yshift_rads);

  double truexcenter_sigma = xcenter_sigma + (shiftx_nn_ * nndist_);
  double reflxcenter_sigma = truexcenter_sigma - xprd*floor(truexcenter_sigma / xprd);

  double trueycenter_sigma = ycenter_sigma + (shifty_nn_ * nndist_);
  double reflycenter_sigma = trueycenter_sigma - yprd*floor(trueycenter_sigma / yprd);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double xdist_to_center = x[i][0] - reflxcenter_sigma;
      double ydist_to_center = x[i][1] - reflycenter_sigma;

      // Look to closest periodic image
      if (xdist_to_center > xprd/2) xdist_to_center -= xprd;
      if (xdist_to_center <-xprd/2) xdist_to_center += xprd;
      if (ydist_to_center > yprd/2) ydist_to_center -= yprd;
      if (ydist_to_center <-yprd/2) ydist_to_center += yprd;


      //choice of corrugation
      if (0) {
        //potential is cropped sines to make cusps
        cx = sin(fmod(xdist_to_center*radspersigma_x/2., M_PI) );
        cy = sin(fmod(ydist_to_center*radspersigma_y/2., M_PI) );
        sx = .5*cos(fmod(xdist_to_center*radspersigma_x/2., M_PI) );
        sy = .5*cos(fmod(ydist_to_center*radspersigma_y/2., M_PI) );
      } else {
        // potential is cosines.  force is negative derivative (sines).
        // origin is the envelope center
        cx= cos(xdist_to_center*radspersigma_x);
        cy= cos(ydist_to_center*radspersigma_y);
        sx= sin(xdist_to_center*radspersigma_x);
        sy= sin(ydist_to_center*radspersigma_y);
      }

      //if ((x[i][0] < 2) && (x[i][1] < 2)) 
      //  printf("x %f y %f cx cy sx sy %f %f %f %f\n", x[i][0], x[i][1], cx, cy, sx, sy);
      rad_sq = (pow(xdist_to_center, 2.0) + 
                pow(ydist_to_center, 2.0));

      // Choice of enevelope
      if (0) {
      // Envelope is a Hertzian profile
        double p0Hertz=2.0 / M_PI * sqrt(aradsq) / r_overEstar; // fixed march 6 
	double norm_to_latpot_coupling = 1; // assumed, but is probably nonlinear
        if (rad_sq < aradsq) envelope = p0Hertz * sqrt(1.0 - rad_sq/aradsq);
        else envelope = 0;
      } else {
      // Envelope is disk with height=1 and cosine roll-off over 1/2 nn spacing
        if (rad_sq < aradsq) envelope = 1;
        else if (rad_sq < arad_border_sq) 
  	  envelope = 0.5 + 0.5*cos(M_PI*(sqrt(rad_sq)-acontact_sigma)/border_dist);
        else envelope = 0;
      }
      // Envelope is diamond-tiled disk with height=1 so V=0 at edge... not implemented

      // foriginal[0] = "potential energy" for added force
      // foriginal[123] = force on atoms before extra force added
      //if (logfile) if (i == 500) fprintf(logfile," adf:i500prev:f[i][0]=%f\n", f[i][0]);
      foriginal[0] += cx*envelope/radspersigma_x + cy*envelope/radspersigma_y; // added feb 12 2014
      foriginal[1] += f[i][0];
      foriginal[2] += f[i][1];
      f[i][0] += sx*envelope;
      f[i][1] += sy*envelope;
      //if (logfile) if (i == 500) fprintf(logfile,"envelope %f at atom x%f y%f and force (%f %f)\n", envelope, x[i][0], x[i][1], sx*envelope, sy*envelope);

      //if ((envelope != 0) && (x[i][0] < 333.))  printf("x %f y %f\n", x[i][0], x[i][1]); 
       // {printf("cx cy sx sy %f %f %f %f\n", cx, cy, sx, sy);}

      //egrav += cx*envelope + cy*envelope;

    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamic::min_post_force(int vflag)
{
  //if (logfile) fprintf(logfile, "ADF calling post force from min_post_force\n");
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */
 // There is a slight error in this calc due to the roll-off at the patch edge
double FixAddForceDynamic::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) { // eflag is set to zero each time post_force called
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForceDynamic::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
