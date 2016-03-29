#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pointers.h"
#ifndef NO_LAMMPS
#include "update.h"
#endif

#include "gfmd_solver_dynamic.h"

/* ----------------------------------------------------------------------
 * Default dynamic implementation using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

GFMDSolverDynamic::GFMDSolverDynamic(LAMMPS *lmp, int narg, int *iarg,
                                     char **arg)
  : GFMDSolverFFT(lmp)
{
  strcpy(name, "dynamic");

  dump_ = NULL;

  if (*iarg > narg-4) {
    error->all(FLERR,"solver default/dynamic: Expected number of layers, layer "
               "spacing, mass and dissipation arguments.");
  }

  char *endptr;

  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected nglob argument.");
  nglob_ = strtol(arg[*iarg], &endptr, 10);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to "
            "integer.", arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected delta argument.");
  delta_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected mass argument.");
  mass_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected gamma argument.");
  gamma_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;

  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "dump")) {
      dump_ = fopen("gfmd_solver_dynamic.out", "w");
      (*iarg)++;
    }
  }

  if (screen)
    fprintf(screen, "Using %i dynamic elastic layers spaced at %f.\n"
	    "Mass of each atom within these layers is %f.\n"
	    "Dynamics is damped with a damping constant of %f*|q|.\n",
	    nglob_, delta_, mass_, gamma_);
  if (logfile)
    fprintf(logfile, "Using %i dynamic elastic layers spaced at %f.\n"
	    "Mass of each atom within these layers is %f.\n"
	    "Dynamics is damped with a damping constant of %f*|q|.\n",
	    nglob_, delta_, mass_, gamma_);

  n_ = NULL;
  u_ = NULL;
  v_ = NULL;
  f_ = NULL;

  q_ = NULL;
  dyn_U0_ = NULL;
  dyn_U_ = NULL;
  dyn_V_ = NULL;
}


GFMDSolverDynamic::~GFMDSolverDynamic()
{
  if (dump_)
    fclose(dump_);

  if (n_) {
    for (int i = 0; i < nxy_loc; i++) {
      memory->destroy(u_[i]);
      memory->destroy(v_[i]);
      memory->destroy(f_[i]);
    }

    destroy_complex_buffer(u0_);
    destroy_complex_buffer(f0_);

    delete n_;
    delete u_;
    delete v_;
    delete f_;
  }

  if (q_) {
    memory->destroy(q_);

    memory->destroy(dyn_U0_);
    memory->destroy(dyn_U_);
    memory->destroy(dyn_V_);
  }
}


void GFMDSolverDynamic::set_grid_size(int in_nx, int in_ny, int in_ndof)
{
  GFMDSolverFFT::set_grid_size(in_nx, in_ny, in_ndof);

  n_ = new int[nxy_loc];
  u_ = new double_complex**[nxy_loc];
  v_ = new double_complex**[nxy_loc];
  f_ = new double_complex**[nxy_loc];

  u0_ = create_complex_buffer("GFMDSolverDynamic::u0");
  f0_ = create_complex_buffer("GFMDSolverDynamic::f0");

  for (int i = 0; i < nxy_loc; i++) {
    n_[i] = nglob_;

    memory->create(u_[i], n_[i], ndof, "GFMDSolverDynamic::u");
    memory->create(v_[i], n_[i], ndof, "GFMDSolverDynamic::v");
    memory->create(f_[i], n_[i], ndof, "GFMDSolverDynamic::f");

    memset(u_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
    memset(v_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
    memset(f_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
  }
}


void GFMDSolverDynamic::set_kernel(StiffnessKernel *kernel, bool normalize)
{
  memory->create(q_, nxy_loc, "GFMDSolverDynamic::q");

  memory->create(dyn_U0_, nxy_loc, ndof_sq, "GFMDSolverDynamic::U0");
  memory->create(dyn_U_, nxy_loc, ndof_sq, "GFMDSolverDynamic::U");
  memory->create(dyn_V_, nxy_loc, ndof_sq, "GFMDSolverDynamic::V");
  memory->create(phi, nxy_loc, ndof_sq, "GFMDSolverDynamic::phi");

  double inv_nxny = 1./double(nx*ny);

  int idq = 0;
  for (int i = xlo_loc; i <= xhi_loc; i++) {
    double qx = (i <= int((nx)/2)) ?
      (2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
    for (int j = ylo_loc; j <= yhi_loc; j++) {
      double qy = (j <= int((ny)/2)) ? 
        (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

      q_[idq] = sqrt( qx*qx + qy*qy );

      kernel->get_dynamical_matrices(qx, qy, dyn_U0_[idq], dyn_U_[idq],
                                     dyn_V_[idq]);

      // this is where our GFunctions^-1 are imported from
      kernel->get_stiffness_matrix(qx, qy, phi[idq]);

      if (i == 0 && j == 0) {
        enforce_phi0_sum_rule(ndof, phi[idq], 1e-15);
      }

      if (normalize) {
        // divide phi_q by (nx*ny) to reduce float operation after FFT
        for (int idim=0; idim<ndof_sq; idim++) {
          // premptively divide for later convenience
          dyn_U0_[idq][idim] *= inv_nxny;
          dyn_U_[idq][idim] *= inv_nxny;
          dyn_V_[idq][idim] *= inv_nxny;
          phi[idq][idim] *= inv_nxny;
        }
      }

      idq++;
    }
  }
}


double GFMDSolverDynamic::post_force(void *input_buffer_ptr,
				     void *output_buffer_ptr,
				     char *dump_prefix)
{
  double **input_buffer = static_cast<double**>(input_buffer_ptr);
  double **output_buffer = static_cast<double**>(output_buffer_ptr);

  /*
   * u(r) -> u(q)
   */
  fft_forward(input_buffer, u0_);

  /*
   * Dump information to file
   */
  if (dump_prefix)
    dump(dump_prefix, u0_);

  verlet_step1();
  double epot = energy_and_forces();
  double ekin = verlet_step2();

  if (dump_)
    fprintf(dump_, "%e %e  %e %e\n", epot0_, epot1_, ekin0_, ekin1_);

  /*
   * f(q) -> f(r)
   */ 
  fft_reverse(f0_, output_buffer);

  return epot + ekin;
}


double GFMDSolverDynamic::energy_and_forces()
{
  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
  double epot0 = 0.0, epot1 = 0.0;
  memset(f0_[0], 0, nxy_loc*ndof*sizeof(double_complex));
  for (int idq = 0; idq < nxy_loc; idq++) {
    int ni = n_[idq];
    double epot = 0.0;
    double_complex *ui0 = u0_[idq];
    double_complex *fi0 = f0_[idq];
    double_complex **ui = u_[idq];
    double_complex **fi = f_[idq];

    double_complex *U0 = dyn_U0_[idq];
    double_complex *U  = dyn_U_[idq];
    double_complex *V  = dyn_V_[idq];

    //    printmat(ndof, V);
    //    printf("ui = %f %f  %f %f  %f %f\n", ui0[0], ui0[1], ui0[2]);

    /* clear all forces */
    memset(fi[0], 0, ni*ndof*sizeof(double_complex));

    /* top layer */
    mat_mul_sub_vec(ndof, U0, ui0,   fi0);
    mat_mul_sub_vec(ndof, V,  ui[0], fi0);
    epot -= creal(conj_vec_dot_vec(ndof, ui0, fi0));

    /* intermediate layers */
    conj_mat_mul_sub_vec(ndof, V, ui0,   fi[0]);
    //    printf("fi = %f %f  %f %f  %f %f\n", fi[0][0], fi[0][1], fi[0][2]);
    mat_mul_sub_vec     (ndof, U, ui[0], fi[0]);
    mat_mul_sub_vec     (ndof, V, ui[1], fi[0]);
    epot -= creal(conj_vec_dot_vec(ndof, ui[0], fi[0]));
    for (int idn = 1; idn < ni-1; idn++) {
      conj_mat_mul_sub_vec(ndof, V, ui[idn-1], fi[idn]);
      mat_mul_sub_vec     (ndof, U, ui[idn],   fi[idn]);
      mat_mul_sub_vec     (ndof, V, ui[idn+1], fi[idn]);
      epot -= creal(conj_vec_dot_vec(ndof, ui[idn], fi[idn]));
    }

    /* bottom layer */
    conj_mat_mul_sub_vec(ndof, V,        ui[ni-2], fi[ni-1]);
    mat_mul_sub_vec     (ndof, phi[idq], ui[ni-1], fi[ni-1]);
    epot -= creal(conj_vec_dot_vec(ndof, ui[ni-1], fi[ni-1]));

    //    printf("%e %e %e\n", ui[3][0], ui[3][1], ui[3][2]);

    /* damping term */
    double gamma_qi = gamma_*q_[idq];
    double_complex **vi = v_[idq];
    for (int idn = 0; idn < ni-1; idn++) {
      for (int j = 0; j < ndof; j++) {
        fi[idn][j] -= gamma_qi*vi[idn][j];
      }
    }

    if (idq == 0) 
      epot0 += epot;
    else
      epot1 += epot;
  }

  epot0_ = epot0;
  epot1_ = epot1;

  return 0.5*(epot0+epot1);
}


void GFMDSolverDynamic::verlet_step1()
{
#ifndef NO_LAMMPS
  double dt = update->dt;
#else
  double dt = 0.1;
#endif
  double d2t = dt*dt;

  for (int idq = 0; idq < nxy_loc; idq++) {
    double_complex **ui = u_[idq];
    double_complex **vi = v_[idq];
    double_complex **fi = f_[idq];

    /* Integration of idn = 0 is treated explicitly in the MD simulation */
    for (int idn = 0; idn < n_[idq]; idn++) {
      for (int idof = 0; idof < ndof; idof++) {

        ui[idn][idof] += vi[idn][idof] * dt + 0.5 * fi[idn][idof] / mass_ * d2t;

        vi[idn][idof] +=                      0.5 * fi[idn][idof] / mass_ * dt;
      }
    }
  }
}


double GFMDSolverDynamic::verlet_step2()
{
#ifndef NO_LAMMPS
  double dt = update->dt;
#else
  double dt = 0.1;
#endif

  double ekin0 = 0.0, ekin1 = 0.0;
  for (int idq = 0; idq < nxy_loc; idq++) {
    double ekin = 0.0;
    double_complex **vi = v_[idq];
    double_complex **fi = f_[idq];

    /* Integration of idn = 0 is treated explicitly in the MD simulation */
    for (int idn = 0; idn < n_[idq]; idn++) {
      for (int idof = 0; idof < ndof; idof++) {
        vi[idn][idof] += 0.5 * fi[idn][idof] / mass_ * dt;
        ekin += 0.5*mass_*creal(conj(vi[idn][idof])*vi[idn][idof]);
      }
    }

    if (idq == 0)
      ekin0 = ekin;
    else
      ekin1 = ekin;
  }

  ekin0_ = ekin0;
  ekin1_ = ekin1;

  return ekin0+ekin1;
}


double GFMDSolverDynamic::memory_usage()
{
  double bytes = 0.0;

  // u, v, f
  for (int i = 0; i < nxy_loc; i++) {
    bytes += 3*n_[i]*ndof*sizeof(double_complex);
  }
  // dyn_U0, dyn_U, dyn_V
  bytes += 3*nxy_loc*ndof_sq*sizeof(double_complex);

  return bytes + GFMDSolverFFT::memory_usage();
}
