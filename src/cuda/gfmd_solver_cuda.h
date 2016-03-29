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
