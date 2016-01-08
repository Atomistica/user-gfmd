#ifndef GFMD_SOLVER_FFT_H
#define GFMD_SOLVER_FFT_H

#include "gfmd_solver.h"
#include "linearalgebra.h"
#include "pointers.h"

namespace LAMMPS_NS {

class GFMDSolverFFT : public GFMDSolver {
 public:
  GFMDSolverFFT(LAMMPS *);
  virtual ~GFMDSolverFFT();

  virtual void set_grid_size(int, int, int);

  virtual double memory_usage();

  virtual void dump_stiffness();
  virtual void dump_greens_function();

  virtual void fft_forward(double **, double_complex **, 
			   double_complex **fac=NULL);
  virtual void fft_reverse(double_complex **, double **,
			   double_complex **fac=NULL);

  double_complex **phi;

 protected:
  class FFT3d *fft;
  double *fft_data;

  void dump(char *, double_complex **);
};

}

#endif
