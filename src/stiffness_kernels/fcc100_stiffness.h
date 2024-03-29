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

StiffnessKernelStyle(fcc100,FCC100StiffnessKernel)

#else

#ifndef __FCC100_STIFFNESS_H
#define __FCC100_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class FCC100StiffnessKernel : public StiffnessKernel {
 public:
    FCC100StiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  virtual ~FCC100StiffnessKernel();

  virtual void set_parameters(double *);

  virtual void get_per_layer_dynamical_matrices(double, double, 
						double_complex **,
						double *, double *);
  virtual void get_force_at_gamma_point(double *);

  virtual void dump_info(FILE *);

 private:
  /*
   * Surface normal
   */
  int z_;

  /*
   * Additional phase multiplier (for debugging purposes)
   */
  double mulphase_;

  /*
   * Forces added to each surface layer
   */
  double *linf_;

  /*
   * Spring constants for nearest- and next-nearest neighbors
   * k is d2V/dr^2_ij at equilibrium r_ij  (Eq A4)
   * kappa is dV/dr_ij * 1/r_ij at equilibrium r_ij (Eq A6)
   */
  double k_[3], kappa_[3];
};

 
 class FCC100x2StiffnessKernel : public StiffnessKernel {
  public:
    FCC100x2StiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
   ~FCC100x2StiffnessKernel();
 
   void get_per_layer_dynamical_matrices(double, double, double_complex **,
 					double *, double *);
 
  private:
   /*
    * Spring constants for nearest- and next-nearest- and third-nearest-neighbors
    */
   double k_[3];
 };

}

#endif

#endif
