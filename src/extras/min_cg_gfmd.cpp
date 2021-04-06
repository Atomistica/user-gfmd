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

#include "math.h"
#include "string.h"
#include "mpi.h"
#include "min_cg_gfmd.h"
#include "atom.h"
#include "fix_gfmd.h"
#include "fix_minimize.h"
#include "modify.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXATOMS 0x7FFFFFFF

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinCGGFMD::MinCGGFMD(LAMMPS *lmp) : MinLineSearch(lmp)
{
  gfmd = NULL;
}

/* ---------------------------------------------------------------------- */

void MinCGGFMD::setup_style()
{
  MinLineSearch::setup_style();

  int i = modify->find_fix("gfmd");
  if (i < 0)
    error->all(FLERR,"minimize cg/gfmd: Could not find gfmd fix.");

  gfmd = (FixGFMD*) modify->fix[i];

  fix_minimize->add_vector(3);
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinCGGFMD::reset_vectors()
{
  MinLineSearch::reset_vectors();

  z = fix_minimize->request_vector(3);
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

int MinCGGFMD::iterate(int maxiter)
{
  int i,fail,ntimestep;
  double beta,gg,dot[3],dotall[3];

  if (!gfmd)
    error->all(FLERR,"minimize cg/gfmd: gfmd fix not set.");

  // we only work with atomic degrees of freedom

  if (nextra_atom > 0 || nextra_global > 0)
    error->all(FLERR,"min cg/gfmd: Cannot work with degrees of freedom other "
	       "than atomic (implement it!).");

  // initialize working vectors

  for (i = 0; i < nvec; i++) h[i] = g[i] = fvec[i];

  // precondition gradient

  for (i = 0; i < nvec; i++)
    z[i] = 0.0;
  gfmd->prec_gradient(fvec, z);
  neval++;

  gg = fnorm_sqr();

  for (int iter = 0; iter < maxiter; iter++) {
    ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    fail = (this->*linemin)(ecurrent,alpha_final);    
    if (fail) return fail;

    // precondition gradient

    for (i = 0; i < nvec; i++)
      z[i] = 0.0;
    gfmd->prec_gradient(fvec, z);
    neval++;

    // Now: fvec is gradient
    //      g is old gradient
    //      z is preconditioned gradient

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) < 
	update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    dot[0] = dot[1] = dot[2] = 0.0;
    for (i = 0; i < nvec; i++) {
      dot[0] += fvec[i]*fvec[i];
      dot[1] += z[i]*fvec[i];
      dot[2] += z[i]*g[i];
    }

    MPI_Allreduce(dot, dotall, 3, MPI_DOUBLE, MPI_SUM, world);

    if (dotall[0] < update->ftol*update->ftol) return FTOL;

    // update new search direction h from new f = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be Fletcher-Reeves
    // reinitialize CG every ndof iterations by setting beta = 0.0

    // fvec is current gradient
    // g is old gradient
    // h is search direction

    beta = MAX(0.0,(dotall[1] - dotall[2])/gg);
    gg = dotall[1];

    for (i = 0; i < nvec; i++) {
      g[i] = fvec[i];
      h[i] = z[i] + beta*h[i];
    }

    // reinitialize CG if new search direction h is not downhill

    dot[0] = 0.0;
    for (i = 0; i < nvec; i++) dot[0] += g[i]*h[i];
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);

    if (dotall[0] <= 0.0) {
      for (i = 0; i < nvec; i++) h[i] = z[i];
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
