/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

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
#include <numeric>

#include <string.h>
#include <stdlib.h>

#include "gfmd_solver_static.h"

#include "pointers.h"

#include "comm.h"
#include "domain.h"
#include "linearalgebra.h"
#include "memory.h"
#include "mpi.h"
#include "fft3d_wrap.h"

#include "gfmd_misc.h"

using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

/* ----------------------------------------------------------------------
 * Default static implementation using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

GFMDSolverStatic::GFMDSolverStatic(LAMMPS *lmp, int narg, int *iarg, char **arg)
  : GFMDSolverFFT(lmp)
{
  strcpy(name, "static");

  linf_ = NULL;
  dump_ = false;
  Gn0_ = NULL;

  if (narg > 0 && *iarg < narg) {
    if (!strcmp(arg[*iarg], "dump")) {
      dump_ = true;
      (*iarg)++;
    }
  }
}


GFMDSolverStatic::~GFMDSolverStatic()
{
  destroy_complex_buffer(q_buffer_);
  if (phi)
    destroy_complex_operator_buffer(phi);
  if (linf_)
    delete [] linf_;
  if (Gn0_)
    memory->destroy(Gn0_);
}


void GFMDSolverStatic::set_grid_size(int in_nx, int in_ny, int in_ndof)
{
  GFMDSolverFFT::set_grid_size(in_nx, in_ny, in_ndof);

  q_buffer_ = create_complex_buffer("GFMDSolverStatic::q_buffer");
  phi = NULL;

  linf_ = new double[nu_];
}


void GFMDSolverStatic::set_kernel(StiffnessKernel *kernel, bool normalize)
{
  if (screen && comm->me == 0)
    fprintf(screen, "USER-GFMD: Computing stiffness matrices...\n");

  /*
   * Get stiffness matrix
   */
  if (phi)
    destroy_complex_operator_buffer(phi);
  phi = create_complex_operator_buffer("GFMDSolverStatic::phi");
  fill_phi_buffer(ndof, nx, xlo_loc, xhi_loc, ny, ylo_loc, yhi_loc,
                  kernel, phi, normalize, error);

  /*
   * Get linear force contributions
   */
  kernel->get_force_at_gamma_point(linf_);

  /*
   * Sanity check: Do forces sum to zero?
   */
  if (std::abs(std::accumulate(linf_, linf_+nu_, 0.0)) > 1e-6) {
      fprintf(screen, "USER-GFMD: Warning: Forces do not sum to zero at the "
              "surface. Sum is = %e. Individual values are",
              std::accumulate(linf_, linf_+nu_, 0.0));
      fprintf(screen, " %e", linf_[0]);
      for (int i = 1; i < nu_; i++) {
          fprintf(screen, ", %e", linf_[i]);
      }
      fprintf(screen, ".\n");
  }

  /*
   * If dump is true we need to compute all Gn0
   */
  if (dump_) {
    height_ = kernel->get_height();
    memory->create(Gn0_, height_, ndof_sq, "GFMDSolverStatic::Gn0");
    for (int i = 0; i < height_; i++) {
      kernel->get_Gn0(i, 0.0, 0.0, Gn0_[i]);
    }
  }

  if (screen && comm->me == 0)
    fprintf(screen, "USER-GFMD: ...done\n");
}


/*
 * FIXME! This dump_prefix stuff is not pretty, but it was the most
 * straighforward way to do that without the need to redo the FFT.
 * Maybe the below function should be split into a fft_forward, apply, and
 * fft_reverse part.
 */
double GFMDSolverStatic::post_force(void *input_buffer_ptr,
				    void *output_buffer_ptr,
				    char *dump_prefix)
{
  double epot;

  double **input_buffer = static_cast<double**>(input_buffer_ptr);
  double **output_buffer = static_cast<double**>(output_buffer_ptr);

  /*
   * u(r) -> u(q)
   */
  fft_forward(input_buffer, q_buffer_);

  /*
   * Store q=0 component
   */
#ifdef HAVE_C99
  double u0_loc[ndof];
  double_complex F_q[ndof];
#else
  double *u0_loc = new double[ndof];
  double_complex *F_q = new double_complex[ndof];
#endif
  memset(u0_loc, 0, ndof*sizeof(double));
  if (xlo_loc == 0 && ylo_loc == 0) {
    //double fac = sqrt(nx*ny);
    double fac = 1.0;
    for (int idim = 0; idim < ndof; idim++) {
      u0_loc[idim] = fac*creal(q_buffer_[0][idim]);
    }
  }
  MPI_Allreduce(u0_loc, u0, ndof, MPI_DOUBLE, MPI_SUM, world);

  /*
   * Dump information to file
   */
  if (dump_prefix)
    dump(dump_prefix, q_buffer_);

  /*
   * Potential energy
   */
  epot = 0.0;

  /*
   * Compute gamma point energy
   */
  if (gammai_ >= 0) {
    for (int i = 0; i < nu_; i++) {
      epot -= 2*linf_[i]*creal(q_buffer_[gammai_][3*i+2]);
    }
  }

  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
  for (int idq = 0; idq < nxy_loc; idq++) {
    mat_mul_vec(ndof, phi[idq], q_buffer_[idq], F_q);
    for (int idim = 0; idim < ndof; idim++) {
      epot += creal(conj(F_q[idim]) * q_buffer_[idq][idim]);
      q_buffer_[idq][idim] = -F_q[idim];
    }
  }  

  /*
   * Add force to gamma point
   */
  if (gammai_ >= 0) {
    if (dump_) {
      FILE *f = fopen("u0.out", "w");
      for (int i = 0; i < height_; i++) {
        mat_mul_vec(ndof, Gn0_[i], q_buffer_[gammai_], F_q);
        fprintf(f, "%i", i);
        for (int j = 0; j < ndof; j++) {
          fprintf(f,"  %e %e", creal(F_q[j]), cimag(F_q[j]));
        }
        fprintf(f, "\n");
      }
      fclose(f);
    }

    for (int i = 0; i < nu_; i++) {
      q_buffer_[gammai_][3*i + 2] += linf_[i];
    }
  }

  /*
   * E = 1/2 U D U - f0*U0
   */
  epot *= 0.5;

  /*
   * f(q) -> f(r)
   */
  fft_reverse(q_buffer_, output_buffer);

#ifndef HAVE_C99
  delete [] u0_loc;
  delete [] F_q;
#endif

  return epot;
}



void GFMDSolverStatic::prec_gradient(double *cavg, double **g, double **gP)
{
  /*
   * g(r) -> g(q)
   */
  fft_forward(g, q_buffer_);

  /*
   * Apply preconditioner
   */
  precondition_gradient<DEF_G>(ndof, nxy_loc, cavg, phi, q_buffer_, error);

  /*
   * gP(q) -> gP(r)
   */ 
  fft_reverse(q_buffer_, gP);
}



