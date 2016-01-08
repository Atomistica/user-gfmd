#ifndef __FCC100_STIFFNESS_H
#define __FCC100_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

using namespace std;

namespace LAMMPS_NS {

class FCC100StiffnessKernel : public StiffnessKernel {
 public:
    FCC100StiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~FCC100StiffnessKernel();

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
