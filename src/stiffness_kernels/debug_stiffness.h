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
