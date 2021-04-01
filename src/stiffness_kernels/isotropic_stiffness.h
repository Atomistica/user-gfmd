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

#ifdef STIFFNESS_KERNEL_CLASS

StiffnessKernelStyle(isotropic,IsotropicStiffnessKernel)
StiffnessKernelStyle(isotropic/z,IsotropicStiffnessKernel)

#else

#ifndef __ISOTROPIC_STIFFNESS_H
#define __ISOTROPIC_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class IsotropicStiffnessKernel : public StiffnessKernel {
 public:
  IsotropicStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  virtual ~IsotropicStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *,
                                    double_complex dU=0.0);

 private:
  double nu_, myu_, gamma_;
};


class IsotropicZStiffnessKernel : public StiffnessKernel {
 public:
  IsotropicZStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  virtual ~IsotropicZStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *); 

 private:
  double nu_, myu_, gamma_;
};

}

#endif

#endif
