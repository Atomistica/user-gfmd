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
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "netcdfcpp.h"

#include "force.h"
#include "pointers.h"
#ifndef NO_LAMMPS
#include "update.h"
#endif

#include "gfmd_analyzer.h"

/* ----------------------------------------------------------------------
 * Default dynamic implementation using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

GFMDAnalyzer::GFMDAnalyzer(LAMMPS *lmp, int narg, int *iarg, char **arg)
  : GFMDSolverFFT(lmp)
{
  strcpy(name, "analyzer");

  dump_ = NULL;
  dump_frame_ = 0;

  if (*iarg > narg-1) {
    error->all(FLERR,"analyzer: Expected number of layers arguments.");
  }

  char *endptr;

  if (*iarg >= narg)
    error->all(FLERR, "analyzer: Expected number of layers argument.");
  n_ = strtol(arg[*iarg], &endptr, 10);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "analyzer: Could not convert '%s' to integer.", arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;

  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "dump")) {
      dump_ = new NcFile("gfmd_analyzer.nc", NcFile::Replace);
      (*iarg)++;
      if (*iarg < narg) {
	dump_interval_ = strtod(arg[*iarg], &endptr);
	if (endptr == arg[*iarg]) {
	  char errstr[1024];
	  sprintf(errstr, "analyzer: Could not convert '%s' to "
		  "double.", arg[*iarg]);
	  error->all(FLERR,errstr);
	}
	dump_counter_ = dump_interval_;
	(*iarg)++;
      }
      else {
	error->all(FLERR,"analyzer: Expected dump interval.");
      }
    }
  }

  u0_ = NULL;
  u1_ = NULL;
  u2_ = NULL;
  v_ = NULL;

  q_ = NULL;
  qx_ = NULL;
  qy_ = NULL;
  linf_ = NULL;
  dyn_U0_ = NULL;
  dyn_U_ = NULL;
  dyn_V_ = NULL;

  ekin_per_q_ = NULL;
  epot_per_q_ = NULL;

  ekin_per_layer_ = NULL;
  epot_per_layer_ = NULL;
}


GFMDAnalyzer::~GFMDAnalyzer()
{
  if (dump_) {
    delete dump_;
  }

  if (q_) {
    memory->destroy(u0_);
    memory->destroy(u1_);
    memory->destroy(u2_);
    memory->destroy(v_);

    memory->destroy(q_);
    memory->destroy(qx_);
    memory->destroy(qy_);

    delete [] linf_;

    memory->destroy(dyn_U0_);
    memory->destroy(dyn_U_);
    memory->destroy(dyn_V_);
  }

  if (ekin_per_q_) {
    memory->destroy(ekin_per_q_);
    memory->destroy(epot_per_q_);

    memory->destroy(ekin_per_layer_);
    memory->destroy(epot_per_layer_);
  }
}


void GFMDAnalyzer::set_grid_size(int in_nx, int in_ny, int in_ndof)
{
  GFMDSolverFFT::set_grid_size(in_nx, in_ny, in_ndof);

  nxy_ = nx*ny;
  inv_nxy_ = 1.0/(nx*ny);

  memory->create(u0_, MAX(1, nxy_loc), ndof, "GFMDAnalyzer::u0");
  memory->create(u1_, MAX(1, nxy_loc), ndof, "GFMDAnalyzer::u0");
  memory->create(u2_, MAX(1, nxy_loc), ndof, "GFMDAnalyzer::u0");
  memory->create(v_, MAX(1, nxy_loc), ndof, "GFMDAnalyzer::u0");

  memory->create(ekin_per_q_, nxy_loc, "GFMDAnalyzer::ekin");
  memory->create(epot_per_q_, nxy_loc, "GFMDAnalyzer::epot");

  memory->create(ekin_per_layer_, n_, "GFMDAnalyzer::ekin");
  memory->create(epot_per_layer_, n_, "GFMDAnalyzer::epot");

  /*
   * Prepare output file structure
   */
  if (dump_) {
    nc_xdim_ = dump_->add_dim("x", nx);
    nc_ydim_ = dump_->add_dim("y", ny);
    nc_ndim_ = dump_->add_dim("layer", n_);
    nc_framedim_ = dump_->add_dim("frame");

    nc_epot_per_q_ = dump_->add_var("epot_per_q", ncDouble, nc_framedim_,
				    nc_xdim_, nc_ydim_);
    nc_epot_per_layer_ = dump_->add_var("epot_per_layer", ncDouble, 
					nc_framedim_, nc_ndim_);
    nc_ekin_per_q_ = dump_->add_var("ekin_per_q", ncDouble, nc_framedim_,
				    nc_xdim_, nc_ydim_);
    nc_ekin_per_layer_ = dump_->add_var("ekin_per_layer", ncDouble,
					nc_framedim_, nc_ndim_);
  }
}


void GFMDAnalyzer::set_kernel(StiffnessKernel *kernel, bool normalize)
{
  kernel_ = kernel;

  memory->create(q_, nxy_loc, "GFMDAnalyzer::q");
  memory->create(qx_, nxy_loc, "GFMDAnalyzer::qx");
  memory->create(qy_, nxy_loc, "GFMDAnalyzer::qy");
  linf_ = new double[nu_];

  memory->create(dyn_U0_, nxy_loc, ndof_sq, "GFMDAnalyzer::U0");
  memory->create(dyn_U_, nxy_loc, ndof_sq, "GFMDAnalyzer::U");
  memory->create(dyn_V_, nxy_loc, ndof_sq, "GFMDAnalyzer::V");

  int idq = 0;
  for (int i = xlo_loc; i <= xhi_loc; i++) {
    double qx = (i <= int((nx)/2)) ?
      (2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
    for (int j = ylo_loc; j <= yhi_loc; j++) {
      double qy = (j <= int((ny)/2)) ? 
        (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

      q_[idq] = sqrt( qx*qx + qy*qy );
      qx_[idq] = qx;
      qy_[idq] = qy;

      kernel->get_dynamical_matrices(qx, qy, dyn_U0_[idq], dyn_U_[idq],
				     dyn_V_[idq]);

      idq++;
    }
  }

  for (int i = 0; i < nu_; i++) {
    kernel->get_force_at_gamma_point(linf_);
  }
}


template<typename T> void rotate_right(T &a, T &b, T &c)
{
  T tmp = c;
  c = b;
  b = a;
  a = tmp;
}


double GFMDAnalyzer::join(double **u_on_grid, double **v_on_grid,
			  char *dump_prefix)
{
  /*
   * Clear all energies
   */
  memset(ekin_per_q_, 0, nxy_loc*sizeof(double));
  memset(epot_per_q_, 0, nxy_loc*sizeof(double));
  memset(ekin_per_layer_, 0, n_*sizeof(double));
  memset(epot_per_layer_, 0, n_*sizeof(double));

  /*
   * u(r) -> u(q), v(r) -> v(q) and analyze
   */
  fft_forward(&u_on_grid[0], u1_);

  for (int idn = 0; idn < n_; idn++) {
    if (idn < n_-1) {
      fft_forward(&u_on_grid[ndof*(idn+1)], u2_);
    }
    fft_forward(&v_on_grid[ndof*idn], v_);
    analyze_layer(u0_, u1_, u2_, v_, idn);
    rotate_right(u2_, u1_, u0_);
  }

  /*
   * clear v_on_grid, because this is supposed to contain the forces on exit
   */
  memset(&v_on_grid[0][0], 0, ndof*n_*nxy_*sizeof(double));

  /*
   * Print per layer energies
   */
#if 0
  for (int i = 0; i < n_; i++)
    printf("%e ", epot_per_layer_[i]);
  printf("\n");
#endif

  /*
   * Sum potential and kinetic energy
   */
  double epot = 0.0;
#if 1
  for (int idq = 0; idq < nxy_loc; idq++) {
    epot += epot_per_q_[idq];
  }
#endif
#if 0
  for (int idn = 0; idn < n_; idn++) {
    epot += epot_per_layer_[idn];
  }
#endif

  if (dump_) {
    dump_counter_--;
    if (dump_counter_ <= 0) {
      nc_ekin_per_q_->set_cur(dump_frame_, 0, 0);
      nc_ekin_per_q_->put(ekin_per_q_, 1, nx, ny);
      nc_epot_per_q_->set_cur(dump_frame_, 0, 0);
      nc_epot_per_q_->put(epot_per_q_, 1, nx, ny);
      nc_ekin_per_layer_->set_cur(dump_frame_, 0);
      nc_ekin_per_layer_->put(ekin_per_layer_, 1, n_);
      nc_epot_per_layer_->set_cur(dump_frame_, 0);
      nc_epot_per_layer_->put(epot_per_layer_, 1, n_);

      dump_counter_ = dump_interval_;
      dump_frame_++;
    }
  }

  return 0.5*epot;
}


void GFMDAnalyzer::analyze_layer(double_complex **u0, double_complex **u1, 
				 double_complex **u2, double_complex **v,
				 int idn)
{
  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q) and compute potential
   * energy 0.5 x F(q) x U(q) per q vector. Do also compute kinetic energy
   * 0.5 x V(q) x V(q) per q vector.
   */
  for (int idq = 0; idq < nxy_loc; idq++) {
    double_complex *ui0 = u0[idq];
    double_complex *ui1 = u1[idq];
    double_complex *ui2 = u2[idq];
    double_complex *vi  = v[idq];

    double_complex *U0  = dyn_U0_[idq];
    double_complex *U   = dyn_U_[idq];
    double_complex *V   = dyn_V_[idq];

    //double_complex fi[ndof_sq];
    double_complex *fi  = v[idq];
   
    /* kinetic energy */
    double ekin = creal(conj_vec_dot_vec(ndof, vi, vi));

#if 0
    if (idq == 0) {
      printf("%i: ui0, ui1, ui2 = %e, %e, %e\n", idn, creal(ui0[2]),
	     creal(ui1[2]), creal(ui2[2]));
    }
#endif

    /* potential energy */
    double epot = 0.0;
    if (idn == 0) {
      /* top layer */
      //mat_mul_sub_vec(ndof, U0, ui1, fi);
      //mat_mul_sub_vec(ndof, V,  ui2, fi);
      mat_mul_vec    (ndof, U0, ui1, fi);
      mat_mul_add_vec(ndof, V,  ui2, fi);
      epot = creal(conj_vec_dot_vec(ndof, ui1, fi));
    }
    else if (idn == n_-1) {
      /* bottom layer */
      //conj_mat_mul_sub_vec(ndof, V,        ui0, fi);
      //conj_mat_mul_sub_vec(ndof, U,        ui1, fi);
      conj_mat_mul_vec    (ndof, V,        ui0, fi);
      conj_mat_mul_add_vec(ndof, U,        ui1, fi);
      epot = creal(conj_vec_dot_vec(ndof, ui1, fi));
    }
    else {
      /* intermediate layers */
      //conj_mat_mul_sub_vec(ndof, V, ui0, fi);
      //mat_mul_sub_vec     (ndof, U, ui1, fi);
      //mat_mul_sub_vec     (ndof, V, ui2, fi);
      conj_mat_mul_vec(ndof, V, ui0, fi);
      mat_mul_add_vec (ndof, U, ui1, fi);
      mat_mul_add_vec (ndof, V, ui2, fi);
      epot = creal(conj_vec_dot_vec(ndof, ui1, fi));
    }

    ekin *= inv_nxy_;
    epot *= inv_nxy_;

    ekin_per_q_[idq] += ekin;
    epot_per_q_[idq] += epot;

    ekin_per_layer_[idn] += ekin;
    epot_per_layer_[idn] += epot;
  }

  /* linear terms */
  if (gammai_ >= 0 && idn == 0) {
    double epot = 0.0;
    for (int iu = 0; iu < nu_; iu++) {
      epot -= 2*linf_[iu]*creal(u1[gammai_][2+3*iu]);
    }

    epot_per_q_[gammai_] += epot;
    epot_per_layer_[idn] += epot;
  }
}


double GFMDAnalyzer::memory_usage()
{
  double bytes = 0.0;

  // u, v
  bytes += 2*nxy_loc*ndof*sizeof(double_complex);
  // dyn_U0, dyn_U, dyn_V
  bytes += 3*nxy_loc*ndof_sq*sizeof(double_complex);
  // ekin_per_q_, epot_per_q_
  bytes += 2*nxy_loc*sizeof(double);
  // ekin_per_layer_, epot_per_layer_
  bytes += 2*n_*sizeof(double);

  return bytes + GFMDSolverFFT::memory_usage();
}
