#ifndef GFMD_SOLVER_STATIC_H
#define GFMD_SOLVER_STATIC_H

#include "gfmd_solver_fft.h"
#include "linearalgebra.h"
#include "pointers.h"

namespace LAMMPS_NS {

class GFMDSolverStatic : public GFMDSolverFFT {
 public:
  GFMDSolverStatic(LAMMPS *, int, int *, char **);
  ~GFMDSolverStatic();

  virtual void set_grid_size(int, int, int);
  virtual void set_kernel(StiffnessKernel *, bool normalize=true);

  virtual double post_force(void *, void *, char *);

  virtual void prec_gradient(double *, double **, double **);

 protected:
  /*
   * Linear force contribution from kernel
   */
  double *linf_;

  /*
   * Dump bulk displacements to file
   */
  bool dump_;

  /*
   * The Green's function Gn0 elements
   */
  int height_;
  double_complex **Gn0_;

  /*
   * Work buffer
   */
  double_complex **q_buffer_;
};

}

#endif
