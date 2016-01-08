#ifndef GFMD_SOLVER_DYNAMIC_H
#define GFMD_SOLVER_DYNAMIC_H

#include "gfmd_solver_fft.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class GFMDSolverDynamic : public GFMDSolverFFT {
 public:
  GFMDSolverDynamic(LAMMPS *, int, int *, char **);
  ~GFMDSolverDynamic();

  void set_grid_size(int, int, int);
  void set_kernel(StiffnessKernel *, bool normalize=true);

  double post_force(void *, void *, char *);

  double memory_usage();

 protected:
  /* value of the reciprocal space vector */
  double *q_;

  /* interaction matrices */
  double_complex **dyn_U0_, **dyn_U_, **dyn_V_;

  /* number of elements in chain per q-vector */
  int *n_, nglob_;

  /* interlayer spacing */
  double delta_;

  /* mass of the atoms */
  double mass_;

  /* damping constants */
  double gamma_;

  /* displacements and forces */
  double_complex **u0_, **f0_;
  double_complex ***u_, ***v_, ***f_;

  /* energies in q=0 and other modes */
  double ekin0_, ekin1_, epot0_, epot1_;

  /* dump file */
  FILE *dump_;

  /* energy_and_forces returns potential, verlet_step 2 return kinetic energy */
  double energy_and_forces();
  void verlet_step1();
  double verlet_step2();
};

}

#endif
