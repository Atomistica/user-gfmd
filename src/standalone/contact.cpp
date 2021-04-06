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
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "contact.h"

#include "fix_contact_map.h"
#include "gfmd_misc.h"
#include "gfmd_stiffness.h"
#include "live_view.h"
#include "nc_traj_io.h"

#include "../REV"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define MAXLEN 1024
#define MAXARGS 64

enum { MODE_NONE, MODE_MIN, MODE_COMPUTE_ONLY };

Contact::Contact(const char *fn) : Pointers()
{
  error = new Error();
  memory = new Memory(error);
  comm = new Comm();
  domain = new Domain();
  screen = stdout;
  logfile = fopen("log.contact", "w");

  strcpy(dump_fn, "traj.nc");

  if (screen)
    fprintf(screen, "# CONTACT r%i, compiled for %i degrees of freedom\n",
	    REV, NDOF);
  if (logfile)
    fprintf(logfile, "# CONTACT r%i, compiled for %i degrees of freedom\n",
	    REV, NDOF);

  restart_fn = NULL;
  outevery = 1;
  delta = -1.0;
  for (int idim = 0; idim < NDOF; idim++)
    fix_center_[idim] = false;

  map = NULL;
  kernel = NULL;
#ifdef HAVE_SDL
  live_view = NULL;
#endif

  mode_ = MODE_NONE;

  nx_ = ny_ = 1;

  rho_tol = 0.01;
  dmax = 0.01;

  ftol = 0.1;
  ntol = -1;
  atol = 0.001;
  max_iter = 500;

  nsteps_ = 1;
  z0_ = -1.0;
  z1_ = 0.0;

  FILE *f = fopen(fn, "r");
  if (!f)
    error->all("Error opening contact.in");

  char line[MAXLEN];
  while (fgets(line, MAXLEN, f)) {
    int i, nargs, lastarg;
    char *args[MAXARGS];

    lastarg = -1;
    nargs = -1;
    args[0] = &line[0];
    i = 0;
    while (line[i]) {
      if (isspace(line[i])) {
	lastarg = nargs;
	line[i] = 0;
      }
      else {
	if (nargs == lastarg) {
	  nargs++;
	  args[nargs] = &line[i];
	}
      }
      i++;
    }
    nargs++;

    if (args[0][0] != '#') {

      if (nargs > 0) {
	if (!strcmp(args[0], "grid")) {
	  if (nargs != 4)
	    error->all("Syntax: grid <nx> <ny> <z0>");
	  nx_ = atoi(args[1]);
	  ny_ = atoi(args[2]);
	  double z0 = atof(args[3]);

	  domain->set_cell(nx_, ny_, 100.0);
	  atom = new Atom(nx_*ny_, memory);

	  double **x = atom->x;
	  double **v = atom->v;
	  double **xeq = atom->xeq;
	  
	  int k = 0;
	  z_ = z0;
	  for (int i = 0; i < nx_; i++) {
	    for (int j = 0; j < ny_; j++) {
	      x[k][0] = i+0.25;
	      x[k][1] = j+0.25;
	      x[k][2] = z0;
	      v[k][0] = v[k][1] = v[k][2] = 0.0;
	      xeq[k][0] = x[k][0];
	      xeq[k][1] = x[k][1];
	      xeq[k][2] = z0;
	      k++;
	    }
	  }
	}
	else if (!strcmp(args[0], "fix")) {
	  if (!strcmp(args[3], "contact/map")) {
	    if (map)
	      error->all("Contact map already defined.");
	    map = new FixContactMap(this, nargs-1, &args[1]);
	    h0_ = map->get_rms_height();
	    avgh_ = map->get_average_height();
	  }
	  else if (!strcmp(args[3], "gfmd") ||
		   !strcmp(args[3], "gfmd/static")) {

	    int iarg = 4;
	    //prefix = strdup(arg[iarg]);
	    iarg++;

	    if ((iarg+1)>nargs) {
	      error->all("fix gfmd: Insufficient command line options. "
			 "Please provide keyword for stiffness expression.");
	    }

	    iarg++;
	    kernel = stiffness_kernel_factory(args[iarg-1], nargs, &iarg, args,
					      domain, memory, error);

	    if (!kernel) {
	      char errstr[1024];
	      sprintf(errstr, "fix gfmd: Unknown analytic stiffness keyword "
		      "'%s' encountered.", args[iarg-1]);
	      error->all(errstr);
	    }

	  
	    // get kernel dimension
	    if (NDOF != kernel->get_dimension()) {
	      error->all("NDOF != kernel->get_dimension. Please recompile.");
	    }

	    while (iarg < nargs) {
	      if (strcmp(args[iarg],"fix_center") == 0) {
		if ((iarg+4)>nargs)
		  error->all("fix gfmd: fix_center: Insufficient command line "
			     "options");
		for (int i = 0; i < 3; i++) {
		  iarg++;
		  fix_center_[i] = atoi(args[iarg]);
		}
		iarg++;
	      }
	      else {
		char str[1024];
		sprintf(str,"fix gfmd: Unknown fix arg: %s", args[iarg]);
		error->all(str);
		iarg++;
	      }
	    }

	  }
	  else {
	    char errstr[1024];
	    sprintf(errstr, "Unknown fix '%s'.", args[3]);
	    error->all(errstr);
	  }
	}
	else if (!strcmp(args[0], "min")) {
	  if (mode_ != MODE_NONE) {
	    error->all("Please use either *min* or *compute_only*.");
	  }
	  mode_ = MODE_MIN;
	  if (!(nargs >= 4 && nargs <= 6)) {
	    error->all("min <ftol|ntol> <tol> <max_iter> [<atol>] [<dmax>] "
		       "[<rho_tol>]\n");
	  }
	  int iarg = 1;
	  if (!strcmp(args[iarg], "ftol")) {
	    iarg++;
	    ftol = atof(args[iarg]);
	    iarg++;
	  }
	  if (!strcmp(args[iarg], "ntol")) {
	    iarg++;
	    ntol = atoi(args[iarg]);
	    iarg++;
	    atol = atof(args[iarg]);
	    iarg++;
	  }
	  else {
	    error->all("min <ftol|ntol> <tol> <max_iter>\n");
	  }
	  if (iarg < nargs) {
	    max_iter = atoi(args[iarg]);
	    iarg++;
	  }
	  if (iarg < nargs) {
	    dmax = atof(args[iarg]);
	    iarg++;
          }
	  if (iarg < nargs) {
	    rho_tol = atof(args[iarg]);
	    iarg++;
	  }
	}
	else if (!strcmp(args[0], "compute_only")) {
	  if (mode_ != MODE_NONE) {
	    error->all("Please use either *min* or *compute_only*.");
	  }
	  mode_ = MODE_COMPUTE_ONLY;
	}
#ifdef HAVE_SDL
	else if (!strcmp(args[0], "live_view")) {
	  if (nargs != 3) {
	    error->all("live_view <sx> <sy>\n");
	  }

	  int sx = atoi(args[1]);
	  int sy = atoi(args[2]);

	  live_view = new LiveView(nx_, ny_, sx, sy, error);
	}
#endif
	else if (!strcmp(args[0], "ramp")) {
	  if (nargs != 4) {
	    error->all("ramp <nsteps> <z0> <dz>\n");
	  }
	  nsteps_ = atoi(args[1]);
	  z0_ = atof(args[2]);
	  z1_ = atof(args[3]);
	}
	else if (!strcmp(args[0], "thermo")) {
	  if (nargs != 2) {
	    error->all("thermo <output frequency>\n");
	  }
	  outevery = atoi(args[1]);
	}
	else if (!strcmp(args[0], "dump")) {
	  if (nargs != 2) {
	    error->all("dump <trajectory filename>");
	  }
	  strcpy(dump_fn, args[1]);
	}
	else if (!strcmp(args[0], "restart")) {
	  if (nargs != 2) {
	    error->all("restart <trajectory filename>");
	  }
	  restart_fn = strdup(args[1]);
	}
	else {
	  char errstr[1024];
	  sprintf(errstr, "Unknown keyword '%s'.", args[0]);
	  error->all(errstr);
	}
      }
      else {
	error->all("Need at least two keywords per line.");
      }

    }
  }
  fclose(f);

#ifdef _OPENMP
  int nthreads = 1;
  nthreads = omp_get_max_threads();

  if (screen)
    fprintf(screen, "# OpenMP: Using %i threads\n", nthreads);
  if (logfile)
    fprintf(logfile, "# OpenMP: Using %i threads\n", nthreads);

  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#endif

  nxy_ = nx_*ny_;
  nq_ = nx_*(ny_/2+1);

  u_r_ = new double*[NDOF];
  f_r_ = new double*[NDOF];

  u_q_ = new double_complex*[NDOF];
  f_q_ = new double_complex*[NDOF];

  fft_forward_u = new fftw_plan[NDOF];
  fft_reverse_u = new fftw_plan[NDOF];
  fft_forward_f = new fftw_plan[NDOF];
  fft_reverse_f = new fftw_plan[NDOF];

  fftq_ = (double_complex*) fftw_malloc(sizeof(double_complex) * nq_);
  fftr_ = (double*) fftq_;

  printf("FFTW3 setup...\n");
  unsigned flags = FFTW_MEASURE;
  for (int idim = 0; idim < NDOF; idim++) {
    u_r_[idim] = (double*)
      fftw_malloc(sizeof(double) * nxy_);
    f_r_[idim] = (double*)
      fftw_malloc(sizeof(double) * nxy_);

    u_q_[idim] = (double_complex*)
      fftw_malloc(sizeof(double_complex) * nq_);
    f_q_[idim] = (double_complex*)
      fftw_malloc(sizeof(double_complex) * nq_);

    /* Note that r2c and c2r overwrites the input */
    fft_forward_u[idim] =
      fftw_plan_dft_r2c_2d(nx_, ny_,
			   fftr_,
			   reinterpret_cast<fftw_complex*>(u_q_[idim]),
			   flags);
    fft_reverse_u[idim] = 
      fftw_plan_dft_c2r_2d(nx_, ny_,
			   reinterpret_cast<fftw_complex*>(fftq_),
			   u_r_[idim],
			   flags);
    fft_forward_f[idim] = 
      fftw_plan_dft_r2c_2d(nx_, ny_,
			   fftr_,
			   reinterpret_cast<fftw_complex*>(f_q_[idim]),
			   flags);
    fft_reverse_f[idim] = 
      fftw_plan_dft_c2r_2d(nx_, ny_,
			   reinterpret_cast<fftw_complex*>(fftq_),
			   f_r_[idim],
			   flags);
  }
  printf("...done\n");

  memory->create(phi_, nq_, NDOF_SQ, "Contact::phi");
  fill_phi_buffer(NDOF, nx_, 0, nx_-1, ny_, 0, ny_/2, kernel, phi_, false,
		  error);
}


Contact::~Contact()
{
  delete kernel;
  delete map;
#ifdef HAVE_SDL
  if (live_view)
    delete live_view;
#endif

  for (int idim = 0; idim < NDOF; idim++) {
    fftw_destroy_plan(fft_forward_u[idim]);
    fftw_destroy_plan(fft_reverse_u[idim]);
    fftw_destroy_plan(fft_forward_f[idim]);
    fftw_destroy_plan(fft_reverse_f[idim]);

    fftw_free(u_r_[idim]);
    fftw_free(f_r_[idim]);

    fftw_free(u_q_[idim]);
    fftw_free(f_q_[idim]);
  }

  fftw_free(fftq_);

  delete fft_forward_u;
  delete fft_reverse_u;
  delete fft_forward_f;
  delete fft_reverse_f;

  delete u_r_;
  delete f_r_;

  delete u_q_;
  delete f_q_;

  delete error;
  delete memory;
  if (atom)
    delete atom;
  delete comm;
  delete domain;
  fclose(logfile);

  if (restart_fn)
    free(restart_fn);

#ifdef _OPENMP
  fftw_cleanup_threads();
#endif
}


void Contact::ur_to_uq()
{
  for (int idim = 0; idim < NDOF; idim++) {
    for (int ixy = 0; ixy < nxy_; ixy++) {
      fftr_[ixy] = u_r_[idim][ixy];
    }
    fftw_execute(fft_forward_u[idim]);
  }
}

void Contact::uq_to_ur()
{
  for (int idim = 0; idim < NDOF; idim++) {
    for (int iq = 0; iq < nq_; iq++) {
      fftq_[iq] = u_q_[idim][iq];
    }
    fftw_execute(fft_reverse_u[idim]);
  }
}

void Contact::fr_to_fq()
{
  for (int idim = 0; idim < NDOF; idim++) {
    for (int ixy = 0; ixy < nxy_; ixy++) {
      fftr_[ixy] = f_r_[idim][ixy];
    }
    fftw_execute(fft_forward_f[idim]);
  }
}

void Contact::fq_to_fr()
{
  for (int idim = 0; idim < NDOF; idim++) {
    for (int iq = 0; iq < nq_; iq++) {
      fftq_[iq] = f_q_[idim][iq];
    }
    fftw_execute(fft_reverse_f[idim]);
  }
}


void Contact::xr_to_uq()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **xeq = atom->xeq;

#pragma omp parallel for			\
  default(none)					\
  firstprivate(n)				\
  shared(x, xeq)
  for (int k = 0; k < n; k++) {
    if (NDOF == 1) {
      u_r_[0][k] = x[k][2] - xeq[k][2];
    }
    else {
      for (int idim = 0; idim < NDOF; idim++) {
	u_r_[idim][k] = x[k][idim] - xeq[k][idim];
      }
    }
  }

  ur_to_uq();
}


void Contact::uq_to_xr()
{
  double **x = atom->x;
  double **f = atom->f;
  double **xeq = atom->xeq;

  /*
   * Do reverse FFT and convert displacements into positions
   */
  uq_to_ur();
  double uz = 0.0;
  if (NDOF == 1) {
    for (int i = 0; i < nxy_; i++) {
      x[i][2] = xeq[i][2] + u_r_[0][i]/nxy_;
      f[i][2] = 0.0;
    }
  }
  else {
    for (int i = 0; i < nxy_; i++) {
      uz += u_r_[2][i]/nxy_;
      for (int idim = 0; idim < NDOF; idim++) {
	x[i][idim] = xeq[i][idim] + u_r_[idim][i]/nxy_;
	f[i][idim] = 0.0;
      }
    }
  }
}


void Contact::wall_interaction()
{
  double **f = atom->f;

  /*
   * Convert displacements to positions
   */
  uq_to_xr();

  /*
   * Compute wall forces
   */
  map->pre_force(0);

  /*
   * Transpose force array and do forward FFT
   */
  frep = 0.0;
  fatt = 0.0;
  if (NDOF == 1) {
    for (int i = 0; i < nxy_; i++) {
      double fz = f[i][2];
      f_r_[0][i] = fz;
      if (fz < 0.0) {
	frep -= fz;
      }
      else {
	fatt += fz;
      }
    }
  }
  else {
    for (int i = 0; i < nxy_; i++) {
      for (int idim = 0; idim < NDOF; idim++) {
	f_r_[idim][i] = f[i][idim];
      }

      double fz = f[i][NDOF-1];
      if (fz < 0.0) {
	frep -= fz;
      }
      else {
	fatt += fz;
      }
    }
  }

  fr_to_fq();
}


double Contact::energy_and_forces(bool fixavg[NDOF], double external_force)
{
  wall_interaction();

  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
  epot_wall_ = map->compute_scalar();
  if (NDOF == 1) {
    curv_wall[0] = fabs(map->curvt[2])/nxy_;   // zz
  }
  else {
    curv_wall[0] = fabs(map->curvt[0])/nxy_;   // xx
    curv_wall[4] = fabs(map->curvt[1])/nxy_;   // yy 
    curv_wall[8] = fabs(map->curvt[2])/nxy_;   // zz
    curv_wall[1] = fabs(map->curvt[3])/nxy_;   // xy
    curv_wall[3] = fabs(map->curvt[3])/nxy_;   // yx
    curv_wall[2] = fabs(map->curvt[4])/nxy_;   // xz
    curv_wall[6] = fabs(map->curvt[4])/nxy_;   // zx
    curv_wall[5] = fabs(map->curvt[5])/nxy_;   // yz 
    curv_wall[7] = fabs(map->curvt[5])/nxy_;   // zy
  }
  nrep = (int) map->compute_vector(4);
  natt = (int) map->compute_vector(5);

  double epot_loc = 0.0, epot0_loc = 0.0;
  int ny2 = ny_/2;
#pragma omp parallel for			\
  default(none)					\
  firstprivate(ny2)				\
  reduction(+:epot_loc) reduction(+:epot0_loc)
  for (int ix = 0; ix < nx_; ix++) {
    int iq = ix*(ny2+1);
    for (int iy = 0; iy <= ny2; iy++) {

      double_complex ui_q[NDOF], fi_q[NDOF];
      for (int idim = 0; idim < NDOF; idim++) {
	ui_q[idim] = u_q_[idim][iq];
      }

      mat_mul_vec(NDOF, phi_[iq], ui_q, fi_q);

      double e = 0;
      for (int idim = 0; idim < NDOF; idim++) {
	e += creal(conj(fi_q[idim]) * ui_q[idim]);
	f_q_[idim][iq] -= fi_q[idim];
      }

      if (iy == 0) {
	if (ix == 0)
	  epot0_loc += e;
	epot_loc += e;
      }
      else {
	epot_loc += 2*e;
      }

      iq++;

    }
  }

  epot0_el_ = 0.5*epot0_loc/nxy_;
  epot_el_ = 0.5*epot_loc/nxy_;

  for (int idim = 0; idim < NDOF; idim++) {
    if (fixavg[idim]) {
      f_q_[idim][0] = 0.0;
    }
  }

  f_q_[NDOF-1][0] += external_force*nxy_;
  epot_ext_ = -( z_*nxy_+creal(u_q_[NDOF-1][0]) )*external_force;

  return epot_el_+epot_wall_+epot_ext_;
}


double Contact::get_fmaxq()
{
  double fmaxq_sq = 0.0;
  for (int idim = 0; idim < NDOF; idim++) {
    fmaxq_sq += creal(f_q_[idim][0]*conj(f_q_[idim][0]));
  }
  for (int iq = 1; iq < nq_; iq++) {
    double f_sq = 0.0;
    for (int idim = 0; idim < NDOF; idim++) {
      f_sq += creal(f_q_[idim][iq]*conj(f_q_[idim][iq]));
    }
    fmaxq_sq = MAX(fmaxq_sq, f_sq);
  }
  return sqrt(fmaxq_sq);
}


double Contact::get_fmaxr()
{
  fq_to_fr();
  double fmaxr_sq = 0.0;
  for (int idim = 0; idim < NDOF; idim++) {
    fmaxr_sq += (f_r_[idim][0]/nxy_)*(f_r_[idim][0]/nxy_);
  }
  for (int ixy = 1; ixy < nxy_; ixy++) {
    double f_sq = 0.0;
    for (int idim = 0; idim < NDOF; idim++) {
      f_sq += (f_r_[idim][ixy]/nxy_)*(f_r_[idim][ixy]/nxy_);
    }
    fmaxr_sq = MAX(fmaxr_sq, f_sq);
  }
  return sqrt(fmaxr_sq);
}


double Contact::get_fnorm()
{
  double fnorm_sq = 0.0;
  int ny2 = ny_/2;
#pragma omp parallel for			\
  default(none)					\
  firstprivate(ny2)				\
  reduction(+:fnorm_sq)
  for (int ix = 0; ix < nx_; ix++) {
    for (int idim = 0; idim < NDOF; idim++) {
      int iq = ix*(ny2+1);
      for (int iy = 0; iy <= ny2; iy++) {
	double f_sq = creal(f_q_[idim][iq]*conj(f_q_[idim][iq]));

	if (iy == 0) {
	  fnorm_sq += f_sq;
	}
	else {
	  fnorm_sq += 2*f_sq;
	}

	iq++;
      }
    }
  }
  return sqrt(fnorm_sq)/nxy_;
}


void Contact::print_status(char s, int iter, double eprevious, double de,
			   double fmax, double delta, double a, int nrep,
			   int natt, double u0, double f0)
{
  bool print = false;

  if (outevery == 1) {
    print = true;
  }
  else if (s == 'F' || s == 'N') {
    print = true;
  }
  else if (s == 'a') {
    outcounter--;
    if (outcounter <= 0) {
      print = true;
      outcounter = outevery;
    }
  }

  if (print) {
    if (screen) {
      fprintf(screen,
	      "%c %5i   %12.5e  ( %12.5e )   %12.5e    %12.5e %12.5e   "
	      "%8i (%4.3f) %8i (%4.3f) %8i (%4.3f)   %12.5e\n",
	      s, iter, eprevious, de,
	      fmax, delta, a,
	      nrep, float(nrep)/(nxy_),
	      natt, float(natt)/(nxy_),
	      nrep+natt, float(nrep+natt)/(nxy_),
	      z_+u0);
      fflush(screen);
    }

    if (logfile) {
      fprintf(logfile,
	      "%c %5i   %12.5e  ( %12.5e )   %12.5e    %12.5e %12.5e   "
	      "%8i (%4.3f) %8i (%4.3f) %8i (%4.3f)   %12.5e\n",
	      s, iter, eprevious, de,
	      fmax, delta, a,
	      nrep, float(nrep)/(nxy_),
	      natt, float(natt)/(nxy_),
	      nrep+natt, float(nrep+natt)/(nxy_),
	      z_+u0);
      fflush(logfile);
    }

#ifdef HAVE_SDL
    if (live_view)
      live_view->update(f_r_[NDOF-1]);
#endif
  }
}


void Contact::compute(bool fixavg[NDOF], double external_force)
{
  xr_to_uq();
  wall_interaction();
  xr_to_uq();

  double eprevious = energy_and_forces(fixavg, external_force);
}


int Contact::min_tr(bool fixavg[NDOF], double external_force)
{
  int ny2 = ny_/2;
  double nq_sq = nq_*((double) nq_);
  double nxy_sq = nxy_*((double) nxy_);

  int ncounter = ntol;
  int a1tol;

  double_complex **h, **g;

  outcounter = outevery;

  if (screen)
    fprintf(screen, "# max_iter = %i\n", max_iter);
  if (logfile)
    fprintf(logfile, "# max_iter = %i\n", max_iter);
  if (ntol > 0) {
    a1tol = (int) (atol*nxy_/ntol)+1;
    if (screen) {
      fprintf(screen, "# Converging contact area.\n");
      fprintf(screen, "# ntol = %i\n", ntol);
      fprintf(screen, "# atol = %e\n", atol);
      fprintf(screen, "# a1tol = %i\n", a1tol);
    }
    if (logfile) {
      fprintf(logfile, "# Converging contact area.\n");
      fprintf(logfile, "# ntol = %i\n", ntol);
      fprintf(logfile, "# atol = %e\n", atol);
      fprintf(logfile, "# a1tol = %i\n", a1tol);
    }
  }
  else {
    if (screen) {
      fprintf(screen, "# Converging force.\n");
      fprintf(screen, "# ftol = %e\n", ftol);
    }
    if (logfile) {
      fprintf(logfile, "# Converging force.\n");
      fprintf(logfile, "# ftol = %e\n", ftol);
    }
  }

  if (screen) {
    fprintf(screen, "# dmax = %e\n", dmax);
    fprintf(screen, "# rho_tol = %e\n", rho_tol);
  }
  if (logfile) {
    fprintf(logfile, "# dmax = %e\n", dmax);
    fprintf(logfile, "# rho_tol = %e\n", rho_tol);
  }

  memory->create(g, NDOF, nq_, "Contact::h");
  memory->create(h, NDOF, nq_, "Contact::g");

  // let contact/map displace atoms initially

  xr_to_uq();
  wall_interaction();
  xr_to_uq();

  // initialize working vectors

  nrepprevious = -ntol-1;
  nattprevious = -ntol-1;
  double eprevious = energy_and_forces(fixavg, external_force);

  // now fvec contains the negative gradient

  double fmax = get_fnorm();

  // precondition gradient, h is preconditioned gradient

#pragma omp parallel for			\
  default(none)					\
  shared(h)
  for (int iq = 0; iq < nq_; iq++) {
    for (int idim = 0; idim < NDOF; idim++) {
      h[idim][iq] = f_q_[idim][iq];
    }
  }
  precondition_gradient<TR_G>(NDOF, nq_, curv_wall, phi_, h, error);

  // gg = g.gP = g.P^-1.g = gP.P.gP

  double gg = 0.0;
#pragma omp parallel for			\
  default(none)					\
  firstprivate(ny2)				\
  shared(g, h)					\
  reduction(+:gg)
  for (int ix = 0; ix < nx_; ix++) {
    for (int idim = 0; idim < NDOF; idim++) {
      int iq = ix*(ny2+1);
      for (int iy = 0; iy <= ny2; iy++) {
	double m = 2;
	if (iy == 0)  m = 1;
	g[idim][iq] = f_q_[idim][iq];
	gg += creal(m*g[idim][iq]*conj(h[idim][iq]));
	iq++;
      }
    }
  }
  gg /= nxy_sq;

  // initialize trust region radius

  //  if (delta < 0.0) {
    delta = MIN(dmax, sqrt(gg));
    //  }

  for (int iter = 1; iter <= max_iter; iter++) {
    // update trust region and search direction

    double a = MIN(1, delta/sqrt(gg));

    // update positions

#pragma omp parallel for \
  default(none)		 \
  firstprivate(a)	 \
  shared(h)
    for (int iq = 0; iq < nq_; iq++) {
      for (int idim = 0; idim < NDOF; idim++) {
	u_q_[idim][iq] += a*h[idim][iq];
      }
    }

    // evaluate energies and forces

    double ecurrent = energy_and_forces(fixavg, external_force);

    // now: u_q, f_q are current position and negative gradient,
    // g is old negative gradient

    // force tolerance criterion
    
    fmax = get_fnorm();

    if (ntol > 0) {
      if (nrep+natt > 0) {
	if (abs(nrep+natt-nrepprevious-nattprevious) < a1tol &&
	    abs(nrep-nrepprevious) < a1tol &&
	    abs(natt-nattprevious) < a1tol &&
	    fmax < ftol) {
	  ncounter--;
	  if (ncounter <= 0) {
	    print_status('N', iter, eprevious, ecurrent-eprevious, fmax,
			 delta, a, nrep, natt, creal(u_q_[NDOF-1][0])/nxy_,
			 creal(f_q_[NDOF-1][0])/nxy_);

	    memory->destroy(g);
	    memory->destroy(h);

	    return iter;
	  }
	}
	else {
	  ncounter = ntol;
	}
      }
      else {
	ncounter = ntol;
      }
    }
    else {
      if (fmax < ftol) {
	print_status('F', iter, eprevious, ecurrent-eprevious, fmax, delta, a,
		     nrep, natt, creal(u_q_[NDOF-1][0])/nxy_,
		     creal(f_q_[NDOF-1][0])/nxy_);

	memory->destroy(g);
	memory->destroy(h);

	return iter;
      }
    }

    // check whether new state should be accepted

    double rho = (eprevious-ecurrent)/((a-a*a/2)*gg);

    if (rho < rho_tol) {
      // state declined

      delta = a*sqrt(gg)/2;

      // rollback to old position change and old gradient

#pragma omp parallel for			\
  default(none)					\
  firstprivate(a)				\
  shared(g, h)
      for (int iq = 0; iq < nq_; iq++) {
	for (int idim = 0; idim < NDOF; idim++) {
	  u_q_[idim][iq] -= a*h[idim][iq];
	  f_q_[idim][iq] = g[idim][iq];
	}
      }

      print_status('r', iter, eprevious, ecurrent-eprevious, fmax, delta, a,
		   nrep, natt, creal(u_q_[NDOF-1][0])/nxy_,
		   creal(f_q_[NDOF-1][0])/nxy_);
    }
    else {
      // state accepted, update preconditioner

#pragma omp parallel for			\
  default(none)					\
  shared(h)
      for (int iq = 0; iq < nq_; iq++) {
	for (int idim = 0; idim < NDOF; idim++) {
	  h[idim][iq] = f_q_[idim][iq];
	}
      }
      precondition_gradient<TR_G>(NDOF, nq_, curv_wall, phi_, h, error);

      // update delta

      double factor = MIN(0.5, (1.0-delta/(2*sqrt(gg))));
      factor = 2*fabs(ecurrent-eprevious)/(factor*sqrt(gg));
      // This is Christophs original formulation
      //delta = MIN(4*delta, MAX(delta/4, factor));
      delta = MIN(2*delta, MAX(delta/2, factor));

      // store gradient in g, compute new gg = g.gP

      gg = 0.0;
#pragma omp parallel for			\
  default(none)					\
  firstprivate(ny2)				\
  shared(g, h)					\
  reduction(+:gg)
      for (int ix = 0; ix < nx_; ix++) {
	for (int idim = 0; idim < NDOF; idim++) {
	  int iq = ix*(ny2+1);
	  for (int iy = 0; iy <= ny2; iy++) {
	    double m = 2;
	    if (iy == 0)  m = 1;
	    g[idim][iq] = f_q_[idim][iq];
	    gg += creal(m*g[idim][iq]*conj(h[idim][iq]));
	    iq++;
	  }
	}
      }
      gg /= nxy_sq;

      // print status

      print_status('a', iter, eprevious, ecurrent-eprevious, fmax, delta, a,
		   nrep, natt, creal(u_q_[NDOF-1][0])/nxy_,
		   creal(f_q_[NDOF-1][0])/nxy_);

      // store eprevious

      eprevious = ecurrent;
      nrepprevious = nrep;
      nattprevious = natt;
    }
  }

  return 0;
}


void Contact::loop()
{
  double **x = atom->x;
  double **xeq = atom->xeq;

  FILE *f = fopen("contact_loop.out", "w");
  fprintf(f, "# 1:zeq 2:z(avg) 3:N 4:A(rep) 5:A(att) 6:A(total) 7:N(rep) "
	  "8:N(att) 9:N(total) 10:n(rep) 11:n(att) 12:n(total) 13:fmax "
	  "14:niter\n");

  if (screen) {
    fprintf(screen, "# nsteps = %i\n", nsteps_);
    fprintf(screen, "# z0 = %f\n", z0_);
    fprintf(screen, "# z1 = %f\n", z1_);
    fprintf(screen, "# fix_center = %i %i %i\n",
	    fix_center_[0], fix_center_[1], fix_center_[2]);
  }
  if (logfile) {
    fprintf(logfile, "# nsteps = %i\n", nsteps_);
    fprintf(logfile, "# z0 = %f\n", z0_);
    fprintf(logfile, "# z1 = %f\n", z1_);
    fprintf(logfile, "# fix_center = %i %i %i\n",
	    fix_center_[0], fix_center_[1], fix_center_[2]);
  }

  NCTrajIO *in_nc = NULL;
  if (restart_fn) {
    if (screen)
      fprintf(screen, "# Restarting from file %s\n", restart_fn);
    if (logfile)
      fprintf(logfile, "# Restarting from file %s\n", restart_fn);
    in_nc = new NCTrajIO(this);
    in_nc->open_for_read(restart_fn);
  }

  NCTrajIO *dump_nc = new NCTrajIO(this);
  dump_nc->open_for_write(dump_fn);
  dump_nc->write_header(0.0, 0.0, 0.0, 0.0);
  dump_nc->write_data();

  for (int it = 0; it < nsteps_; it++) {
    if (screen)
      fprintf(screen, "\n=== ITERATION %i of %i ===\n", it+1, nsteps_);
    if (logfile)
      fprintf(logfile, "\n=== ITERATION %i of %i ===\n", it+1, nsteps_);

#if 1
    if (nsteps_ == 1) {
      z_ = avgh_+h0_*z0_;
    }
    else {
      z_ = avgh_+h0_*log(exp(z0_)+it*(exp(z1_)-exp(z0_))/(nsteps_-1));
      //z_ = avgh_+log(exp(h0_*z0_)+it*(exp(h0_*z1_)-exp(h0_*z0_))/(nsteps_-1));
      //z_ = z0_ + it*z1_*nx_;
    }
    for (int i = 0; i < nxy_; i++) {
      xeq[i][2] = z_;
    }

    if (screen)
      fprintf(screen, "# z = %f\n", z_);
    if (logfile)
      fprintf(logfile, "# z = %f\n", z_);
#endif

#if 0
    double external_force = z0_+z1_*it;
    if (screen)
      fprintf(screen, "# external_force = %f\n", external_force);
    if (logfile)
      fprintf(logfile, "# external_force = %f\n", external_force);

    int niter = min_tr(external_force);
#endif

    if (in_nc) {
      if (screen)
	fprintf(screen, "# reading restart configuration\n");
      if (logfile)
	fprintf(logfile, "# reading restart configuration\n");

      in_nc->read_frame(it+1);
    }

    int niter;
    if (mode_ == MODE_MIN) {
      niter = min_tr(fix_center_, 0.0);
    }
    else if (mode_ == MODE_COMPUTE_ONLY) {
      niter = 0;
      compute(fix_center_, 0.0);
    }
    else
      error->all("Internal error: Unknown mode_.");
    double fmax = get_fmaxr();
    //fq_to_fr();
    dump_nc->write_header(epot_wall_+epot_el_, epot_wall_, epot0_el_,
			  epot_el_-epot0_el_);
    dump_nc->write_data();

    fprintf(f, "%e  %e  %e  %e %e %e  %e %e %e  %i %i %i   %e   %i   %e %e "
	    "%e\n",
	    z_, z_+creal(u_q_[NDOF-1][0])/nxy_,
	    -creal(f_q_[NDOF-1][0])/nxy_,
	    ((double) nrep)/nxy_, ((double) natt)/nxy_,
	    ((double) (nrep+natt))/nxy_,
	    frep/nxy_, fatt/nxy_, (frep-fatt)/nxy_,
	    nrep, natt, nrep+natt,
	    fmax, niter,
	    epot_wall_, epot_el_, epot_ext_);
    fflush(f);
  }

  delete dump_nc;

  if (in_nc)
    delete in_nc;

  fclose(f);
}


int main(int argc, char *argv[])
{
  if (argc != 2) {
    printf("Syntax: contact <control file>\n");
    return 999;
  }

  Contact *contact = new Contact(argv[1]);

  contact->loop();

  delete contact;
}
