#ifndef __FCC100FT_STIFFNESS_H
#define __FCC100FT_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

using namespace std;

namespace LAMMPS_NS {

class FCC100FTStiffnessKernel : public StiffnessKernel {
 public:
    FCC100FTStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  ~FCC100FTStiffnessKernel();

  virtual void set_parameters(double *);

  virtual void get_per_layer_dynamical_matrices(double, double, 
						double_complex **,
						double *, double *);
  virtual void get_force_at_gamma_point(double *);

  virtual void dump_info(FILE *);

  int disttoconnecttype(double);
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
   */
  double k_[3], kappa_[3];
};

 
 class FCC100FTx2StiffnessKernel : public StiffnessKernel {
  public:
    FCC100FTx2StiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
   ~FCC100FTx2StiffnessKernel();
 
   void get_per_layer_dynamical_matrices(double, double, double_complex **,
 					double_complex *);
 
  private:
   /*
    * Spring constants for nearest- and next-nearest- and third-nearest-neighbors
    */
   double k_[3];
 };

}

#endif
