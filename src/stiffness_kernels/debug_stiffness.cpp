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
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "debug_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * For debugging purposes only.
 * Assembles synthetic U0, U, V matrices to check if the assembly
 * process work properly.
 * --------------------------------------------------------------------*/

DebugStiffnessKernel::DebugStiffnessKernel(int narg, int *carg, char **arg,
					   Domain *domain, Memory *memory,
					   Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  strcpy(name_, "debug");

  if (*carg >= narg) {
    error_->all(FLERR,"DebugStiffnessKernel::DebugStiffnessKernel: Expected "
		"kernel dimension.");
  }

  nu_ = strtol(arg[*carg], &endptr, 10);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"DebugStiffnessKernel::DebugStiffnessKernel: Error "
		"converting dimension argument to integer.");
  }
  (*carg)++;


  dim_ = ndof_*nu_;
}

DebugStiffnessKernel::~DebugStiffnessKernel()
{
}

void DebugStiffnessKernel::get_per_layer_dynamical_matrices(double qx,
                                                            double qy,
                                                            double_complex **D,
                                                            double *xcshift,
                                                            double *ycshift)
{
  for (int i = 0; i < 2*nu_+1; i++) {
    for (int k = 0; k < ndof_; k++) {
      for (int l = 0; l < ndof_; l++) {
	MEL(ndof_, D[i], k, l) =
	  (i*1000.+k*10.+l*1.) + (i*10000.+k*10.+l*1.)*I;
      }
    }
  }
}


void DebugStiffnessKernel::get_dynamical_matrices(double qx, double qy,
                                                  double_complex *U0,
                                                  double_complex *U,
                                                  double_complex *V,
                                                  double_complex dU)
{
  StiffnessKernel::get_dynamical_matrices(qx, qy, U0, U, V);

  printf("U0: \n");
  printmat(dim_, U0);

  printf("\nU: \n");
  printmat(dim_, U);

  printf("\nV: \n");
  printmat(dim_, U);

  printf("\n");
}


/* ----------------------------------------------------------------------
 * Linear harmonic chain, analytical solution is 1/N where N is the
 * number of layers. For debugging purposes.
 * --------------------------------------------------------------------*/

ChainStiffnessKernel::ChainStiffnessKernel(int narg, int *carg, char **arg,
                                           Domain *domain, Memory *memory,
                                           Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  dim_ = 3;
  strcpy(name_, "chain");

  if (*carg < narg-1) {
    if (!strcmp(arg[*carg], "height")) {
      char *endptr;
      (*carg)++;
      height_ = strtol(arg[*carg], &endptr, 10);
      if (endptr == arg[*carg]) {
	error_->all(FLERR,"ChainStiffnessKernel::ChainStiffnessKernel: Error "
		    "converting height argument to number.");
      }
      (*carg)++;
    }
  }
}

ChainStiffnessKernel::~ChainStiffnessKernel()
{
}

void ChainStiffnessKernel::get_per_layer_dynamical_matrices(double qx,
                                                            double qy,
                                                            double_complex **D,
                                                            double *xcshift,
                                                            double *ycshift)
{
  if (nu_ != 1) {
    error_->all(FLERR,"ChainStiffnessKernel::get_per_layer_dynamical_matrices:"
		"Only works for nu == 1.");
  }

  SquareMatrix<double_complex> U0(dim_, D[0]), U(dim_, D[1]), V(dim_, D[2]);

  U0.fill_with(0.0);
  U.fill_with(0.0);
  V.fill_with(0.0);
 
  U0[0][0] = 1.0;
  U0[1][1] = 1.0;
  U0[2][2] = 1.0;

  U[0][0] = 2.0;
  U[1][1] = 2.0;
  U[2][2] = 2.0;

  V[0][0] = -1.0;
  V[1][1] = -1.0;
  V[2][2] = -1.0;
}

} /* namespace */
