#ifndef __ISOTROPIC_STIFFNESS_H
#define __ISOTROPIC_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

using namespace std;

namespace LAMMPS_NS {

class IsotropicStiffnessKernel : public StiffnessKernel {
 public:
    IsotropicStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~IsotropicStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *,
                                    double_complex dU=0.0);

 private:
  double nu_, myu_, gamma_;
};


class IsotropicZStiffnessKernel : public StiffnessKernel {
 public:
  IsotropicZStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~IsotropicZStiffnessKernel();

  virtual void get_stiffness_matrix(double, double, double_complex *); 

 private:
  double nu_, myu_, gamma_;
};

}

#endif
