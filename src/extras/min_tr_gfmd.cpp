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
#include "mpi.h"
#include "min_tr_gfmd.h"
#include "atom.h"
#include "fix_gfmd.h"
#include "modify.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"

// DEBUG
#include "comm.h"

using namespace LAMMPS_NS;

#define MAXATOMS 0x7FFFFFFF

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinTRGFMD::MinTRGFMD(LAMMPS *lmp) : MinLineSearch(lmp)
{
  gfmd = NULL;
}

/* ---------------------------------------------------------------------- */

void MinTRGFMD::setup_style()
{
  MinLineSearch::setup_style();

  int i = modify->find_fix("gfmd");
  if (i < 0)
    error->all(FLERR,"minimize tr/gfmd: Could not find gfmd fix.");

  gfmd = (FixGFMD*) modify->fix[i];
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

int MinTRGFMD::iterate(int maxiter)
{
  int i,ntimestep;
  double rho_tol = 0.01;

  if (!gfmd)
    error->all(FLERR,"minimize tr/gfmd: gfmd fix not set.");

  // we only work with atomic degrees of freedom

  if (nextra_atom > 0 || nextra_global > 0)
    error->all(FLERR,"min tr/gfmd: Cannot work with degrees of freedom other "
	       "than atomic (implement it!).");

  // initialize working vectors

  eprevious = energy_force(1);
  neval++;
  // now fvec contains the negative gradient

  // precondition gradient, h is preconditioned gradient

  for (i = 0; i < nvec; i++) {
    h[i] = 0.0;
  }

  gfmd->prec_gradient(fvec, h);
  neval++;
  //  printf("%i  %f %f %f  %f %f %f\n", comm->me, fvec[0], fvec[1], fvec[2], h[0], h[1], h[2]);

  // gg = g.gP = g.P^-1.g = gP.P.gP

  double tmp = 0.0, gg;
  for (i = 0; i < nvec; i++) {
    g[i] = fvec[i];
    tmp += g[i]*h[i];
  }
  MPI_Allreduce(&tmp, &gg, 1, MPI_DOUBLE, MPI_SUM, world);

  // initialize trust region radius

  double delta = MIN(dmax, sqrt(gg));

  //  printf("- %i %f %f\n", niter, eprevious, delta);

  for (int iter = 0; iter < maxiter; iter++) {
    ntimestep = ++update->ntimestep;
    niter++;

    // update trust region and search direction

    double a = MIN(1, delta/sqrt(gg));

    // update positions

    for (i = 0; i < nvec; i++) {
      xvec[i] += a*h[i];
    }

    // evaluate energies and forces

    //    printf("A\n");
    ecurrent = energy_force(1);
    neval++;
    //    printf("B\n");

    // now: xvec, fvec are current position and negative gradient,
    // g is old negative gradient

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) < 
	update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion
    
    tmp = 0.0;
    for (i = 0; i < nvec; i++) {
      tmp += fvec[i]*fvec[i];
    }
    double fnorm;
    MPI_Allreduce(&tmp, &fnorm, 1, MPI_DOUBLE, MPI_SUM, world);

    if (fnorm < update->ftol*update->ftol) return FTOL;

    // check whether new state should be accepted

    double rho = (eprevious-ecurrent)/((a-a*a/2)*gg);

    if (rho < rho_tol) {
      // state declined

      delta = a*sqrt(gg)/2;

      // rollback to old position change and old gradient

      for (i = 0; i < nvec; i++) {
	xvec[i] -= a*h[i];
	fvec[i] = g[i];
      }

      //      printf("R %i %f %f %f %f\n", niter, ecurrent, delta, rho, sqrt(fnorm));
    }
    else {
      // state accepted, update preconditioner

      for (i = 0; i < nvec; i++) {
	h[i] = 0.0;
      }

      gfmd->prec_gradient(fvec, h);
      neval++;
      //      printf("%i  %f %f %f  %f %f %f\n", comm->me, fvec[0], fvec[1], fvec[2], h[0], h[1], h[2]);

      // update delta

      double factor = MIN(0.5, (1.0-delta/(2*sqrt(gg))));
      factor = 2*abs(ecurrent-eprevious)/(factor*sqrt(gg));
      //      delta = MIN(4*delta, MAX(delta/4, factor));
      delta = MIN(2*delta, MAX(delta/2, factor));

      // store eprevious

      eprevious = ecurrent;

      // store gradient in g, compute new gg = g.gP

      tmp = 0.0;
      for (i = 0; i < nvec; i++) {
	g[i] = fvec[i];
	tmp += g[i]*h[i];
      }
      MPI_Allreduce(&tmp, &gg, 1, MPI_DOUBLE, MPI_SUM, world);

      // print status

      //      printf("A %i %f %f %f %f %f\n", niter, ecurrent, delta, rho, sqrt(fnorm), gg);
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  return MAXITER;
}
