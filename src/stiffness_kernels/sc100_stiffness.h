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
#ifndef __SC100_STIFFNESS_H
#define __SC100_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

using namespace std;

namespace LAMMPS_NS {

class SC100ExplicitStiffnessKernel : public StiffnessKernel {
 public:
  SC100ExplicitStiffnessKernel(int, int *, char **, Domain *, Memory *,
			       Error *);
  ~SC100ExplicitStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *);
};


class SC100StiffnessKernel : public StiffnessKernel {
 public:
   SC100StiffnessKernel(int, int *, char **,  Domain *, Memory *, Error *);
  ~SC100StiffnessKernel();

  virtual void get_per_layer_dynamical_matrices(double, double, 
						double_complex **,
						double_complex *);
};

}

#endif
