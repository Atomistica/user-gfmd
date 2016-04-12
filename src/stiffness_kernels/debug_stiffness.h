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
#ifndef __DEBUG_STIFFNESS_H
#define __DEBUG_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

using namespace std;

namespace LAMMPS_NS {

class DebugStiffnessKernel : public StiffnessKernel {
 public:
   DebugStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~DebugStiffnessKernel();

  virtual void get_per_layer_dynamical_matrices(double, double,
                                                double_complex **,
                                                double *, double *);
  virtual void get_dynamical_matrices(double, double, double_complex *,
                                      double_complex *, double_complex *,
                                      double_complex dU=0.0);
};


class ChainStiffnessKernel : public StiffnessKernel {
 public:
    ChainStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~ChainStiffnessKernel();

  virtual void get_per_layer_dynamical_matrices(double, double,
                                                double_complex **,
                                                double *, double *);
};

}

#endif
