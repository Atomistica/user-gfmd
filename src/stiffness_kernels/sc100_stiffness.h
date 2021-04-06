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

#ifdef STIFFNESS_KERNEL_CLASS

StiffnessKernelStyle(sc100_explicit,SC100ExplicitStiffnessKernel)
StiffnessKernelStyle(sc100,SC100StiffnessKernel)

#else

#ifndef __SC100_STIFFNESS_H
#define __SC100_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class SC100ExplicitStiffnessKernel : public StiffnessKernel {
 public:
  SC100ExplicitStiffnessKernel(int, int *, char **, Domain *, Memory *,
			       Error *);
  virtual ~SC100ExplicitStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *);
};


class SC100StiffnessKernel : public StiffnessKernel {
 public:
  SC100StiffnessKernel(int, int *, char **,  Domain *, Memory *, Error *);
  virtual ~SC100StiffnessKernel();

  virtual void get_per_layer_dynamical_matrices(double, double, 
						double_complex **,
						double *, double *);
};

}

#endif

#endif
