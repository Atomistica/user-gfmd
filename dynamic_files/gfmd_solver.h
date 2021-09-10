#ifndef GFMD_SOLVER_H
#define GFMD_SOLVER_H

#include "linearalgebra.h"
#include "pointers.h"
#include "surface_stiffness.h"

#define MAX_NDOF 24

namespace LAMMPS_NS {

class GFMDSolver : protected Pointers {
 public:
  GFMDSolver(LAMMPS *);
  virtual ~GFMDSolver();

  virtual void init() { };

  virtual void set_grid_size(int, int, int);
  virtual void set_kernel(StiffnessKernel *, bool normalize=true) = 0;

  /*!
   * The solver can be called asynchronously. Input and output buffers are
   * passed as void *. For on CPU operation they expect a double ** (allocated
   * by the LAMMPS memory object). For GPU operation they expect a double *
   * to the GPU location that holds the memory. FIXME! Harmonize!
   */
  virtual void pre_force(void *, void *) { };
  virtual double post_force(void *, void *, char *) = 0;
  //virtual double post_force(void *, void *, void *, char *) = 0; // pass vxy for dynamic

  // preconditioned gradient
  virtual void prec_gradient(double *, double **, double **) {
    error->all(FLERR,"prec_gradient not implemented for this solver.");
  }

  // preconditioned gradient in reciprocal space
  virtual void prec_gradient_q(double, double **) {
    error->all(FLERR,"prec_gradient_q not implemented for this solver.");
  }

  double **create_double_buffer(const char *);
  void destroy_double_buffer(double **);

  double_complex **create_complex_buffer(const char *);
  void destroy_complex_buffer(double_complex **);

  double_complex **create_complex_operator_buffer(const char *);
  void destroy_complex_operator_buffer(double_complex **);

  int get_xlo_loc() {
    return xlo_loc;
  }
  int get_xhi_loc() {
    return xhi_loc;
  }
  int get_ylo_loc() {
    return ylo_loc;
  }
  int get_yhi_loc() {
    return yhi_loc;
  }
  int get_nxy_loc() {
    return nxy_loc;
  }

  double *get_u0() {
    return u0;
  }

  char *get_name() {
    return name;
  };

  virtual int get_nmin() {
    return 0;
  }
  virtual int get_nmax() {
    return 0;
  }
  virtual int get_lb_flag() {
    return 0;
  }
  
  virtual double memory_usage();

  virtual void dump_stiffness();
  virtual void dump_greens_function();

 protected:
  /* Descriptor for this grid class */
  char name[80];

  /* Parameters */
  int nx, ny, nu_, ndof, ndof_sq;
  double nxtarget;
  double nytarget;
  
  /* Wall indices */
  int *nxwall, *nywall;
  int *fxwall, *fywall;
  
  /* Derived FFT stuff */

  /* my process rank, total number of procs */
  int me, nprocs;

  /* job assignment info */
  int xlo_loc, xhi_loc, nx_loc;
  int ylo_loc, yhi_loc, ny_loc;
  int nxy_loc;

  /* Dynamic with load balancing variables */
  int xlo_loc_ks, xhi_loc_ks, nx_loc_ks;
  int ylo_loc_ks, yhi_loc_ks, ny_loc_ks;
  int nxy_loc_ks, lb_flag_, nmin_, nmax_;
  
  int ntot;
  int *edges_x, *edges_y;

  int nprocs_x, nprocs_y, nprocs_z;
  int nwalls_x, nwalls_y;
  
  /* grid index of the gamma point */
  int gammai_;

  /* q=0 displacement */
  double *u0;

  /* Sum over layers */
  int layer_sum(int,int,int,int,int,int,int,int);
};

GFMDSolver *gfmd_solver_factory(char *, LAMMPS *, int, int *, char **);

}

#endif
