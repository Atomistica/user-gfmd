/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */
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

  virtual double memory_usage();

  virtual void dump_stiffness();
  virtual void dump_greens_function();

 protected:
  /* Descriptor for this grid class */
  char name[80];

  /* Parameters */
  int nx, ny, nu_, ndof, ndof_sq;

  /* Derived FFT stuff */

  /* my process rank, total number of procs */
  int me, nprocs;

  /* job assignment info */
  int xlo_loc, xhi_loc, nx_loc;
  int ylo_loc, yhi_loc, ny_loc;
  int nxy_loc;

  /* grid index of the gamma point */
  int gammai_;

  /* q=0 displacement */
  double *u0;
};

GFMDSolver *gfmd_solver_factory(char *, LAMMPS *, int, int *, char **);

}

#endif
