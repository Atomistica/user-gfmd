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
#ifndef GFMD_SOLVER_CUDA_H
#define GFMD_SOLVER_CUDA_H

#ifdef GFMD_CUFFT
#include <cufft.h>
#endif

#include "gfmd_solver.h"
#include "linearalgebra.h"
#include "pointers.h"

namespace LAMMPS_NS {

class GFMDSolverStaticCUDA : public GFMDSolver {
  public:
    GFMDSolverStaticCUDA(LAMMPS *);
    ~GFMDSolverStaticCUDA();

    virtual void init();

    virtual void set_kernel(StiffnessKernel *, bool normalize=true);

    virtual double post_force(void *, void *, char *);

  protected:
    /*!
     * Linear force contribution from kernel
     */
    double *linf_, *dev_linf_;

    /*!
     * Reciprocal space buffer size (half the size in y-direction)
     */
    int nbuffer_half_loc_;

    /*!
     * CUFFT forward and inverse plans
     */
    cufftHandle fft_plan_forward_, fft_plan_inverse_;

    /*!
     * "Operator" (stiffness matrix) buffer on device
     */
    cufftDoubleComplex *dev_q_operator_;
    size_t dev_q_operator_pitch_;

    /*!
     * Work space on device
     */
    cufftDoubleComplex *dev_q_buffer1_, *dev_q_buffer2_;
    size_t dev_q_buffer_pitch_;

    /*!
     * Accumulator for potential energy
     */
    cufftDoubleReal *dev_epot_buf_;
};

}

#endif
