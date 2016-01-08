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
