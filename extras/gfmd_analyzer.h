#ifndef GFMD_ANALYZER_H
#define GFMD_ANALYZER_H

#include "gfmd_solver_fft.h"
#include "surface_stiffness.h"

#include "netcdfcpp.h"

namespace LAMMPS_NS {

class GFMDAnalyzer : public GFMDSolverFFT {
 public:
  GFMDAnalyzer(LAMMPS *, int, int *, char **);
  ~GFMDAnalyzer();

  void set_grid_size(int, int, int);
  void set_kernel(StiffnessKernel *, bool normalize=true);

  double join(double **, double **, char *);

  double memory_usage();

 protected:
  /* Stiffness kernel */
  StiffnessKernel *kernel_;

  /* number of layers */
  int n_;

  /* value of the reciprocal space vector */
  double *q_, *qx_, *qy_;

  /* linear force terms */
  double *linf_;

  /* interaction matrices */
  double_complex **dyn_U0_, **dyn_U_, **dyn_V_;

  /* number of elements in chain per q-vector */
  bigint nxy_;
  double inv_nxy_;

  /* interlayer spacing */
  double delta_;

  /* mass of the atoms */
  double mass_;

  /* damping constants */
  double gamma_;

  /* displacements and forces */
  double_complex **u0_, **u1_, **u2_, **v_;

  /* energies in each q mode */
  double *ekin_per_q_, *epot_per_q_;

  /* energies in each layer */
  double *ekin_per_layer_, *epot_per_layer_;

  /* dump file, interval and counter */  
  NcFile *dump_;
  int dump_interval_, dump_counter_, dump_frame_;
  NcDim *nc_xdim_, *nc_ydim_, *nc_ndim_, *nc_framedim_;
  NcVar *nc_epot_per_q_, *nc_epot_per_layer_;
  NcVar *nc_ekin_per_q_, *nc_ekin_per_layer_;

  /* energy_and_forces returns potential, verlet_step 2 return kinetic energy */
  void analyze_layer(double_complex **, double_complex **, double_complex **,
		     double_complex **, int);
};

}

#endif
