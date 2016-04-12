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
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <ctype.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#ifndef NO_LAMMPS
#include "update.h"
#include "fft3d_wrap.h"
#endif

#include "fix_wall_map.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { AVERAGE, BOTTOM, TOP };

const char syntax[] =
  "fix wall/map: Illegal fix command, syntax is 'fix ID group-ID "
  "wall/map plane|<mapfile> [bot|top] <z0> <surface normal> "
  "<lj93|lj93/smooth|lj126|exp|expf|dugdale> "
  "[<epsilon> <sigma> <cutoff>] [lambda <short-wavelength cutoff>] "
  "[pendist <distance from surface after penetration>] "
  "[vel <x-velocity> <y-velocity>]'";

#define MODULO(i, n)  ( (i) >= 0 ? (i) % (n) : (n+i) % (n) )

#define round(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))

/* ---------------------------------------------------------------------- */

FixWallMap::FixWallMap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7)
    error->all(FLERR,syntax);

  cutoff_ = -1.0;
  offset_ = 0.0;
  lambda_ = -1.0;
  pendist_ = -1.0;
  firstcall_ = true;

  map_data_ = NULL;
  map_ = NULL;

  x0_ = 0.0;
  y0_ = 0.0;
  vx0_ = 0.0;
  vy0_ = 0.0;

  bool dump_map = false;
  bool lowmem = false;
  bool disable_interpolation = false;
  int iarg = 3;
  char *mapfn = arg[iarg++];

  // read map

  if (!strcmp(mapfn, "flat")) {
    init_flat();
  }
  else if (!strcmp(mapfn, "100")) {
    init_100(iarg, narg, arg);
  }
  else if (!strcmp(mapfn, "111")) {
    init_111(iarg, narg, arg);
  }
  else {
    read_map(mapfn, map_data_);
  }

  // map position

  z0_type_ = AVERAGE;
  if (!strcmp(arg[iarg], "bot")) {
    z0_type_ = BOTTOM;
    iarg++;
  }
  else if (!strcmp(arg[iarg], "top")) {
    z0_type_ = TOP;
    iarg++;
  }

  char *endptr;
  z0_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg]) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Could not convert z0 (= '%s') parameter "
	    "to floating point value. Or did you mean to specify 'bot' or "
	    "'top'?",
            arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;

  // surface normal

  if (!strcmp(arg[iarg], "1")) {
    norm = 1;
  }
  else if (!strcmp(arg[iarg], "-1")) {
    norm = -1;
  }
  else {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Surface normal must be 1 or -1, but "
	    "is '%s'.", arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;

  // type of potential

  if (!strcmp(arg[iarg], "lj93")) {
    iarg++;
    init_lj93(iarg, narg, arg);
  }
  else if (!strcmp(arg[iarg], "lj126")) {
    iarg++;  
    init_lj126(iarg, narg, arg);
  }
  else if (!strcmp(arg[iarg], "lj93/smooth")) {
    iarg++;
    init_lj93smooth(iarg, narg, arg);
  }
  else if (!strcmp(arg[iarg], "exp")) {
    iarg++; 
    init_exp(iarg, narg, arg);
  }
  else if (!strcmp(arg[iarg], "expf")) {
    iarg++;
    init_expf(iarg, narg, arg);
  }
  else if (!strcmp(arg[iarg], "dugdale")) {
    iarg++;
    init_dugdale(iarg, narg, arg);
  }
  else {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Unknown potential type '%s'.",
	    arg[iarg]);
    error->all(FLERR,errstr);
  }

  while (iarg < narg) {
    if (!strcmp(arg[iarg], "lambda")) {
      iarg++;
      if (iarg >= narg)
	error->all(FLERR,"fix contact/map: Expected value after 'lambda' flag.");
      lambda_ = strtod(arg[iarg], &endptr);
      if (endptr == arg[iarg]) {
	char errstr[1024];
	sprintf(errstr, "fix contact/map: Could not convert 'lambda' (= '%s') "
		"parameter to floating point value.", arg[iarg]);
	error->all(FLERR,errstr);
      }
    }
    else if (!strcmp(arg[iarg], "pendist")) {
      iarg++;
      if (iarg >= narg)
	error->all(FLERR,"fix contact/map: Expected value after 'pendist' flag.");
      pendist_ = strtod(arg[iarg], &endptr);
      if (endptr == arg[iarg]) {
	char errstr[1024];
	sprintf(errstr, "fix contact/map: Could not convert 'pendist' (= '%s') "
		"parameter to floating point value.", arg[iarg]);
	error->all(FLERR,errstr);
      }
    }
    else if (!strcmp(arg[iarg], "vel")) {
      iarg++;
      if (iarg >= narg-1)
	error->all(FLERR,"fix contact/map: Expected x- and y-velocity values after "
		   "'vel' flag.");
      vx0_ = strtod(arg[iarg], &endptr);
      if (endptr == arg[iarg]) {
	char errstr[1024];
	sprintf(errstr, "fix contact/map: Could not convert first 'vel' "
		"parameter (= '%s') to floating point value.", arg[iarg]);
	error->all(FLERR,errstr);
      }
      iarg++;
      vy0_ = strtod(arg[iarg], &endptr);
      if (endptr == arg[iarg]) {
	char errstr[1024];
	sprintf(errstr, "fix contact/map: Could not convert second 'vel' "
		"parameter (= '%s') to floating point value.", arg[iarg]);
	error->all(FLERR,errstr);
      }
    }
    else if (!strcmp(arg[iarg], "dump")) {
      dump_map = true;
    }
    else if (!strcmp(arg[iarg], "lowmem")) {
      lowmem = true;
    }
    else if (!strcmp(arg[iarg], "disable_interpolation")) {
      /*
       * Disable spline interpolation, just take value from map
       */
      disable_interpolation = true;
      lowmem = true;
    }
    else {
      char errstr[1024];
      sprintf(errstr, "fix contact/map: Unknown keyword '%s'.", arg[iarg]);
      error->all(FLERR,errstr);
    }
    iarg++;
  }


  if (lambda_ > 0.0) {
    if (!map_data_)
      error->all(FLERR,"Map filtering can only be used with map data file.");
    filter_map(map_data_);
  }
  if (map_data_) {
    displace_map(map_data_);
    if (disable_interpolation) {
      eval_map_func_ = &FixWallMap::eval_raw_data;
    }
    else {
      map_ = new Table2D(nx_loc_, ny_loc_, map_data_, true, lowmem, error,
			 memory);
      eval_map_func_ = &FixWallMap::eval_table2d;
    }
    if (!lowmem) {
      memory->destroy(map_data_);
      map_data_ = NULL;
    }
  }

  // dump map
  if (dump_map) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    int mx = 128;
    int my = 128;
    FILE *f = fopen("map.dump", "w");
    for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
	double v, dvdx, dvdy, dvdxdx, dvdydy, dvdxdy;
	(this->*eval_map_func_)((i+0.5)*xprd/mx, (j+0.5)*yprd/my,
				v, dvdx, dvdy, dvdxdx, dvdydy, dvdxdy);
	fprintf(f, "%20.10e ", v);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }

  // fix behavior

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 6;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  if (fabs(vx0_) > 0.0 || fabs(vy0_) > 0.0) {
    restart_global = 1;
  }

  // other stuff

  force_flag = 0;
  fmap_loc_[0] = fmap_loc_[1] = fmap_loc_[2] = fmap_loc_[3] = fmap_loc_[4] = 
    fmap_loc_[5] = fmap_loc_[6] = 0.0;
  fmap_[0] = fmap_[1] = fmap_[2] = fmap_[3] = fmap_[4] = fmap_[5] = 
    fmap_[6] = 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixWallMap::write_restart(FILE *fp)
{
  int n = 0;
  double list[2];
  list[n++] = x0_;
  list[n++] = y0_;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixWallMap::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  x0_ = static_cast<int> (list[n++]);
  y0_ = static_cast<int> (list[n++]);

  if (screen)
    fprintf(screen, "Read map offset %f %f from restart "
	    "file.\n", x0_, y0_);
  if (logfile)
    fprintf(logfile, "Read map offset %f %f from restart "
	    "file.\n", x0_, y0_);
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_epsilon_and_sigma(int &iarg, int narg, char **arg)
{
  char *endptr;
  if (iarg == narg)  error->all(FLERR,syntax);
  epsilon_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg])
    error->all(FLERR,"fix contact/map: Could not convert first potential "
	       "parameter to floating point value.");
  iarg++;

  if (iarg == narg)  error->all(FLERR,syntax);
  sigma_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg])
    error->all(FLERR,"fix contact/map: Could not convert second potential "
	       "parameter to floating point value.");
  iarg++;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_cutoff(int &iarg, int narg, char **arg)
{
  char *endptr;
  if (iarg == narg)  error->all(FLERR,syntax);
  cutoff_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg])
    error->all(FLERR,"fix contact/map: Could not convert potential cutoff "
	       "parameter to floating point value.");
  iarg++;

  mcutoff_ = 0.01*cutoff_;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_outer_cutoff(int &iarg, int narg, char **arg)
{
  char *endptr;
  inner_cutoff_ = cutoff_;
  if (iarg == narg)  error->all(FLERR,syntax);
  cutoff_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg])
    error->all(FLERR,"fix contact/map: Could not convert second potential cutoff "
	       "parameter to floating point value.");
  iarg++;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_lj93(int &iarg, int narg, char **arg)
{
  init_epsilon_and_sigma(iarg, narg, arg);
  init_cutoff(iarg, narg, arg);

  coeff1_ = 6.0/5.0 * epsilon_ * pow(sigma_,9.0);
  coeff2_ = 3.0 * epsilon_ * pow(sigma_,3.0);
  coeff3_ = 2.0/15.0 * epsilon_ * pow(sigma_,9.0);
  coeff4_ = epsilon_ * pow(sigma_,3.0);
  coeff5_ = 12.0 * epsilon_ * pow(sigma_,9.0);
  coeff6_ = 12.0 * epsilon_ * pow(sigma_,3.0);

  if (cutoff_ > 0.0) {
    double rinv;
    rinv = 1.0/cutoff_;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;

    offset_ = coeff3_*r4inv*r4inv*rinv - coeff4_*r2inv*rinv;
  }

  mc0_ = coeff3_/pow(mcutoff_, 9) - coeff4_/pow(mcutoff_, 3) - offset_;
  mc1_ = - coeff1_/pow(mcutoff_, 10) + coeff2_/pow(mcutoff_, 4);
  mc2_ = coeff5_/pow(mcutoff_, 11) - coeff6_/pow(mcutoff_, 5);

  eval_pot_func_ = &FixWallMap::eval_lj93;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_lj93(double dz, double &e, double &f, double &d2e)
{
  if (dz < mcutoff_) {
    double ddz = dz - mcutoff_;
    e = mc0_ + mc1_*ddz + 0.5*mc2_*ddz*ddz;
    f = - mc1_ - mc2_*ddz;
    d2e = mc2_;
  }
  else {
    double rinv = 1.0/dz;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    double r10inv = r4inv*r4inv*r2inv;
    e = coeff3_*r4inv*r4inv*rinv - coeff4_*r2inv*rinv - offset_;
    f = coeff1_*r10inv - coeff2_*r4inv;
    d2e = coeff5_*r10inv*rinv - coeff6_*r4inv*rinv;  
  }
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_lj126(int &iarg, int narg, char **arg)
{
  init_epsilon_and_sigma(iarg, narg, arg);
  init_cutoff(iarg, narg, arg);

  coeff1_ = 48.0 * epsilon_ * pow(sigma_,12.0);
  coeff2_ = 24.0 * epsilon_ * pow(sigma_,6.0);
  coeff3_ = 4.0 * epsilon_ * pow(sigma_,12.0);
  coeff4_ = 4.0 * epsilon_ * pow(sigma_,6.0);

  if (cutoff_ > 0.0) {
    double r2inv = 1.0/(cutoff_*cutoff_);
    double r6inv = r2inv*r2inv*r2inv;
    offset_ = r6inv*(coeff3_*r6inv - coeff4_);
  }

  eval_pot_func_ = &FixWallMap::eval_lj126;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_lj126(double dz, double &e, double &f, double &d2e)
{
  if (dz < mcutoff_) {
    error->one(FLERR,"fix contact/map: Atom too close: Implement short ranged "
	       "spline for lj126 potential.");
  }

  double rinv = 1.0/dz;
  double r2inv = rinv*rinv;
  double r6inv = r2inv*r2inv*r2inv;
  e = r6inv*(coeff3_*r6inv - coeff4_) - offset_;
  f = r6inv*(coeff1_*r6inv - coeff2_) * rinv;
  d2e = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_lj93smooth(int &iarg, int narg, char **arg)
{
  init_lj93(iarg, narg, arg);
  init_outer_cutoff(iarg, narg, arg);

  double rinv;
  rinv = 1.0/inner_cutoff_;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  double r10inv = r4inv*r4inv*r2inv;

  // make force and first derivative of force continuous
  smcff1_ = coeff1_*r10inv - coeff2_*r4inv;
  smcff2_ = -10*coeff1_*r10inv*rinv + 4*coeff2_*r4inv*rinv;

  double dr = cutoff_ - inner_cutoff_;
  if (dr <= 0.0)
    error->all(FLERR,"fix contact/map: lj93/smooth outer cutoff must be bigger "
	       "than inner cutoff.");
  smcff3_ = -(3*smcff1_ + 2*smcff2_*dr)/(dr*dr);
  smcff4_ = (2*smcff1_ + smcff2_*dr)/(dr*dr*dr);

  smcff0_ = smcff1_*dr + smcff2_*dr*dr/2.0 + smcff3_*dr*dr*dr/3.0 +
    smcff4_*dr*dr*dr*dr/4.0;

  // make energies continuous
  offset_ -= smcff0_;

  eval_pot_func_ = &FixWallMap::eval_lj93smooth;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_lj93smooth(double dz, double &e, double &f,
				    double &d2e)
{
  if (dz < mcutoff_) {
    double ddz = dz - mcutoff_;
    e = mc0_ + mc1_*ddz + 0.5*mc2_*ddz*ddz;
    f = - mc1_ - mc2_*ddz;
    d2e = mc2_;
  }
  else if (dz < inner_cutoff_) {
    double rinv = 1.0/dz;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    double r10inv = r4inv*r4inv*r2inv;
    e = coeff3_*r4inv*r4inv*rinv - coeff4_*r2inv*rinv - offset_;
    f = coeff1_*r10inv - coeff2_*r4inv;
    d2e = coeff5_*r10inv*rinv - coeff6_*r4inv*rinv;
  }
  else {
    double rinv = dz - inner_cutoff_;
    double r2inv = rinv*rinv;
    e = smcff0_ - smcff1_*rinv - smcff2_*r2inv/2.0 -  smcff3_*r2inv*rinv/3.0
      - smcff4_*r2inv*r2inv/4.0;
    f = smcff1_ + smcff2_*rinv + smcff3_*r2inv + smcff4_*r2inv*rinv;
    d2e = -smcff2_ - 2.0*smcff3_*rinv - 3.0*smcff4_*r2inv;
  }
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_exp(int &iarg, int narg, char **arg)
{
  init_epsilon_and_sigma(iarg, narg, arg);
  init_cutoff(iarg, narg, arg);

  coeff1_ = 1.0/sigma_;
  if (cutoff_ > 0.0) {
    offset_ = epsilon_*exp(-cutoff_*coeff1_);
  }

  eval_pot_func_ = &FixWallMap::eval_exp;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_exp(double dz, double &e, double &f, double &d2e)
{
  e  = epsilon_*exp(-dz*coeff1_);
  f  = e*coeff1_;
  d2e = f*coeff1_;
  e -= offset_;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_expf(int &iarg, int narg, char **arg)
{
  init_epsilon_and_sigma(iarg, narg, arg);
  init_cutoff(iarg, narg, arg);

  coeff1_ = 1.0/sigma_;
  if (cutoff_ > 0.0) {
    offset_ = epsilon_*(sigma_+cutoff_);
  }

  eval_pot_func_ = &FixWallMap::eval_expf;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_expf(double dz, double &e, double &f, double &d2e)
{
  d2e = 0.0;
  double tmp = epsilon_*exp(coeff1_*(cutoff_-dz));
  f = tmp - epsilon_;
  d2e = coeff1_*tmp;
  e = sigma_*tmp + epsilon_*dz - offset_;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::init_dugdale(int &iarg, int narg, char **arg)
{
  init_epsilon_and_sigma(iarg, narg, arg);

  cutoff_ = sigma_;
  coeff1_ = epsilon_/cutoff_;

  eval_pot_func_ = &FixWallMap::eval_dugdale;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_dugdale(double dz, double &e, double &f, double &d2e)
{
  d2e = 0.0;
  if (dz > 0.0) {
    e = coeff1_*(cutoff_-dz);
    f = coeff1_;
  }
  else {
    e = 1e12;
    f = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixWallMap::~FixWallMap()
{
  if (map_)
    delete map_;

  if (map_data_)
    memory->destroy(map_data_);
}

/* ---------------------------------------------------------------------- */

int FixWallMap::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= THERMO_ENERGY;
  mask |= MIN_PRE_FORCE;
  mask |= END_OF_STEP;
  nevery = 1;
  return mask;
}

/* ---------------------------------------------------------------------- */

#ifndef NO_LAMMPS
void FixWallMap::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all(FLERR,"fix contact/map: RESPA not yet supported.");
}

/* ---------------------------------------------------------------------- */

void FixWallMap::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    pre_force(vflag);
}
#endif

/* ---------------------------------------------------------------------- */

void FixWallMap::min_setup(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallMap::pre_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (domain->triclinic)
    error->all(FLERR,"fix contact/map: Triclinic boxes not supported.");

  fmap_loc_[0] = fmap_loc_[1] = fmap_loc_[2] = fmap_loc_[3] = fmap_loc_[4] = 
    fmap_loc_[5] = fmap_loc_[6] = 0.0;
  curvt[0] = curvt[1] = curvt[2] = curvt[3] = curvt[4] = curvt[5] = 0.0;

  // integrate offset
  x0_ += vx0_*update->dt;
  y0_ += vy0_*update->dt;

  // for hard wall we also need to reset the positions of ghost atoms,
  // since they won't be communicated after pre_force
  //  if (pot == HARD)
  //    nlocal += atom->nghost;

  double penmax = -1.0;
  bool reset_atom = false;

  double epot = 0.0;
  double fx_sum = 0.0, fy_sum = 0.0, fz_sum = 0.0;
  int ncon = 0, nrep = 0, natt = 0;

  double curvtxx = 0.0, curvtyy = 0.0, curvtzz = 0.0;
  double curvtxy = 0.0, curvtyz = 0.0, curvtzx = 0.0;

#pragma omp parallel for default(none)					\
  shared(mask, f, penmax, x)						\
  firstprivate(nlocal, reset_atom)					\
  reduction(+:epot)							\
  reduction(+:fx_sum) reduction(+:fy_sum) reduction(+:fz_sum)		\
  reduction(+:ncon) reduction(+:nrep) reduction(+:natt)			\
  reduction(+:curvtxx) reduction(+:curvtyy) reduction(+:curvtzz)	\
  reduction(+:curvtxy) reduction(+:curvtyz) reduction(+:curvtzx)
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double z, dz_dx, dz_dy, d2z_dxdx, d2z_dydy, d2z_dxdy;
      (this->*eval_map_func_)(x[i][0], x[i][1], z, dz_dx, dz_dy,
			      d2z_dxdx, d2z_dydy, d2z_dxdy);

      double dz = norm*(x[i][2]-z);

      if (cutoff_ < 0.0 || ( cutoff_ > 0.0 && dz < cutoff_ )) {
	double fz = 0.0, curv = 0.0, epot_loc = 0.0;

	double rinv, r2inv, r4inv, r6inv, r10inv, tmp;

	if (dz < pendist_) {
	  x[i][2] = z+norm*pendist_;
	  dz = pendist_;
	  reset_atom = true;
	}

	(this->*eval_pot_func_)(dz, epot_loc, fz, curv);

        if (fz > 0.0) {
	  nrep++;
        }
        else {
	  natt++;
        }
	ncon++;

	fz *= norm;

	double fx  = -fz*dz_dx;
	double fy  = -fz*dz_dy;

	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	epot += epot_loc;
	fx_sum += fx;
	fy_sum += fy;
	fz_sum += fz;

	curvtxx += fz*d2z_dxdx + curv*dz_dx*dz_dx;    // xx
	curvtyy += fz*d2z_dydy + curv*dz_dy*dz_dy;    // yy
	curvtzz += curv;                              // zz

	curvtxy += fz*d2z_dxdy + curv*dz_dx*dz_dy;    // xy
	curvtyz -= curv*dz_dx;                        // yz
	curvtzx -= curv*dz_dy;                        // zx
      }
    }
  }

  fmap_loc_[0] = epot;
  fmap_loc_[1] = fx_sum;
  fmap_loc_[2] = fy_sum;
  fmap_loc_[3] = fz_sum;
  fmap_loc_[4] = ncon;
  fmap_loc_[5] = nrep;
  fmap_loc_[6] = natt;

  curvt[0] = curvtxx;
  curvt[1] = curvtyy;
  curvt[2] = curvtzz;
  curvt[3] = curvtxy;
  curvt[4] = curvtyz;
  curvt[5] = curvtzx;


  if (penmax > 0.0) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Atom penetrated wall; maximum penetration was %f.", penmax);
    error->one(FLERR,errstr);
  }

  if (reset_atom) {
    if (screen)
      fprintf(screen, "Warning: fix contact/map: Resetting atom position "
	      "because of wall penetration.\n");
    if (logfile)
      fprintf(logfile, "Warning: fix contact/map: Resetting atom position "
	      "because of wall penetration.\n");
  }

  force_flag = 0;
  firstcall_ = false;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ----------------------------------------------------------------------
   return the contact area
------------------------------------------------------------------------- */

double FixWallMap::compute_scalar()
{
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(fmap_loc_, fmap_, 7, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fmap_[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixWallMap::compute_vector(int n)
{
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(fmap_loc_, fmap_, 7, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fmap_[n+1];
}

/* ----------------------------------------------------------------------
   create a surface map that corresponds to a flat plane
------------------------------------------------------------------------- */

void FixWallMap::init_flat()
{
  /*
   * Set evaluation function
   */
  eval_map_func_ = &FixWallMap::eval_flat;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_flat(double x, double y, double &z,
			      double &dz_dx, double &dz_dy,
			      double &d2z_dxdx, double &d2z_dydy,
			      double &d2z_dxdy)
{
  z = z0_;
  dz_dx = 0.0;
  dz_dy = 0.0;
  d2z_dxdx = 0.0;
  d2z_dydy = 0.0;
  d2z_dxdy = 0.0;
}

/* ----------------------------------------------------------------------
   create a surface map that corresponds to the leading Fourier
   components of a (100) surface
------------------------------------------------------------------------- */

// number of grid points per unit cell length
#define RESOLUTION 100
void FixWallMap::init_100(int &iarg, int narg, char **arg)
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;

  /*
   * read lattice constant
   */
  char *endptr;
  double a0 = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg]) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Could not convert lattice constant "
	    "(= '%s') parameter for (111) surface to floating point value.",
	    arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;


  /*
   * read amplitude
   */
  h0_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg]) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Could not convert amplitude (= '%s') "
	    "parameter for (111) surface to floating point value.", arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;

  h0_ /= 2;

  /*
   * Reciprocal lattice vectors
   */
  qx_ = 2*M_PI*round(xprd/a0)/xprd;
  qy_ = 2*M_PI*round(yprd/a0)/xprd;

  if (screen)
    fprintf(screen, "\nCreating (100) surface with lattice parameters "
	    "%f and %f, and %i and %i periods in x- and y-direction, "
	    "respectively.\n", 2*M_PI/qx_, 2*M_PI/qy_,
	    (int) (xprd*qx_/(2*M_PI)), (int) (yprd*qy_/(2*M_PI)));
  if (logfile)
    fprintf(logfile, "\nCreating (100) surface with lattice parameters "
	    "%f and %f, and %i and %i periods in x- and y-direction, "
	    "respectively.\n", 2*M_PI/qx_, 2*M_PI/qy_,
	    (int) (xprd*qx_/(2*M_PI)), (int) (yprd*qy_/(2*M_PI)));

  /*
   * Set evaluation function
   */
  eval_map_func_ = &FixWallMap::eval_100;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_100(double x, double y, double &z,
			     double &dz_dx, double &dz_dy,
			     double &d2z_dxdx, double &d2z_dydy,
			     double &d2z_dxdy)
{
  double dx, dy, cx, cy, sx, sy;

  dx = x - x0_;
  dy = y - y0_;
  cx = cos(qx_*dx);
  cy = cos(qy_*dy);
  sx = sin(qx_*dx);
  sy = sin(qy_*dy);

  z = z0_ + h0_*( cx + cy );
  dz_dx = -h0_*qx_ * sx;
  dz_dy = -h0_*qy_ * sy;
  d2z_dxdx = -h0_*qx_*qx_ * cx;
  d2z_dydy = -h0_*qy_*qy_ * cy;
  d2z_dxdy = 0.0;
}

/* ----------------------------------------------------------------------
   create a surface map that corresponds to the leading Fourier
   components of a (111) surface
------------------------------------------------------------------------- */

void FixWallMap::init_111(int &iarg, int narg, char **arg)
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;

  /*
   * read lattice constant
   */
  char *endptr;
  double a0 = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg]) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Could not convert lattice constant "
	    "(= '%s') parameter for (111) surface to floating point value.",
	    arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;


  /*
   * read amplitude
   */
  h0_ = strtod(arg[iarg], &endptr);
  if (endptr == arg[iarg]) {
    char errstr[1024];
    sprintf(errstr, "fix contact/map: Could not convert amplitude (= '%s') "
	    "parameter for (111) surface to floating point value.", arg[iarg]);
    error->all(FLERR,errstr);
  }
  iarg++;

  h0_ /= 3;

  /*
   * Reciprocal lattice vectors
   */
  double sqrt3 = sqrt(3.0);
  qx_ = 2*M_PI*round(xprd/(a0*sqrt3))/xprd;
  qy_ = 2*M_PI*round(yprd/a0)/yprd;

  if (screen)
    fprintf(screen, "\nCreating (111) surface with lattice parameters "
	    "%f and %f, and %i and %i periods in x- and y-direction, "
	    "respectively.\n", 2*M_PI/(sqrt3*qx_), 2*M_PI/qy_,
	    (int) (xprd*sqrt3*qx_/(2*M_PI)), (int) (yprd*qy_/(2*M_PI)));
  if (logfile)
    fprintf(logfile, "\nCreating (111) surface with lattice parameters "
	    "%f and %f, and %i and %i periods in x- and y-direction, "
	    "respectively.\n", 2*M_PI/(sqrt3*qx_), 2*M_PI/qy_,
	    (int) (xprd*sqrt3*qx_/(2*M_PI)), (int) (yprd*qy_/(2*M_PI)));


  /*
   * Set evaluation function
   */
  eval_map_func_ = &FixWallMap::eval_111;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_111(double x, double y, double &z,
			     double &dz_dx, double &dz_dy,
			     double &d2z_dxdx, double &d2z_dydy,
			     double &d2z_dxdy)
{
  double dx, dy, c1, c2, c3, s1, s2, s3;

  dx = x - x0_;
  dy = y - y0_;
  c1 = cos(qy_*dy);
  c2 = cos(qx_*dx + 0.5*qy_*dy);
  c3 = cos(qx_*dx - 0.5*qy_*dy);
  s1 = sin(qy_*dy);
  s2 = sin(qx_*dx + 0.5*qy_*dy);
  s3 = sin(qx_*dx - 0.5*qy_*dy);

  z = z0_ + h0_*( c1 + c2 + c3 );
  dz_dx = -h0_*qx_*( s2 + s3 );
  dz_dy = -h0_*qy_*( s1 + 0.5*s2 - 0.5*s3 );
  d2z_dxdx = -h0_*qx_*qx_*( c2 + c3 );
  d2z_dydy = -h0_*qy_*qy_*( c1 + 0.25*c2 + 0.25*s3 );
  d2z_dxdy = -0.5*h0_*qx_*qy_*( c2 - c3 );
}

/* ----------------------------------------------------------------------
   allocate map section that resides on the local processor and set 
   dimension variable
------------------------------------------------------------------------- */

#define BUF 16
void FixWallMap::allocate_map(double **&map_data, int in_nx, int in_ny)
{
  nx_ = in_nx;
  ny_ = in_ny;

  /*
   * Check which part of the grid should be stored on this processor
   */
  xlo_loc0_ = (int) (nx_*(domain->sublo[0]-domain->boxlo[0])/domain->xprd);
  xhi_loc0_ = (int) (nx_*(domain->subhi[0]-domain->boxlo[0])/domain->xprd-1);

  ylo_loc0_ = (int) (ny_*(domain->sublo[1]-domain->boxlo[1])/domain->yprd);
  yhi_loc0_ = (int) (ny_*(domain->subhi[1]-domain->boxlo[1])/domain->yprd-1);


  /*
   * Allow some buffer zone at the boundary
   */
  xlo_loc_ = xlo_loc0_ - BUF;
  xhi_loc_ = xhi_loc0_ + BUF;
  ylo_loc_ = ylo_loc0_ - BUF;
  yhi_loc_ = yhi_loc0_ + BUF;


  /*
   * Set exactly to system dimension of only 1 processor
   */
  if (comm->procgrid[0] == 1) {
    xlo_loc_  = 0;
    xhi_loc_  = nx_-1;
    xlo_loc0_ = 0;
    xhi_loc0_ = nx_-1;
  }
  if (comm->procgrid[1] == 1) {
    ylo_loc_  = 0;
    yhi_loc_  = ny_-1;
    ylo_loc0_ = 0;
    yhi_loc0_ = ny_-1;
  }


  /*
   * Local buffer sizes
   */
  nx_loc0_  = xhi_loc0_-xlo_loc0_+1;
  ny_loc0_  = yhi_loc0_-ylo_loc0_+1;

  nx_loc_   = xhi_loc_-xlo_loc_+1;
  ny_loc_   = yhi_loc_-ylo_loc_+1;

  nxy_loc0_ = nx_loc0_*ny_loc0_;
  nxy_loc_  = nx_loc_*ny_loc_;

  memory->create(map_data, nx_loc_, ny_loc_, "FixWallMap::map");
}

/* ----------------------------------------------------------------------
   read a file containing a surface map
------------------------------------------------------------------------- */

void FixWallMap::read_map(char *fn, double **&map_data)
{
  FILE *f = fopen(fn, "r");

  if (!f) {
    char errstr[1024];
    sprintf(errstr, "Could not open contact map file '%s'.", fn);
    error->all(FLERR,errstr);
  }


  /*
   * count number of columns
   */
  nx_ = 0;
  char c = fgetc(f), lastc = ' ';
  while (!feof(f) && c != '\n') {
    if ((isalnum(c) || c == '-' || c == 'e') && lastc == ' ')
      nx_++;
    lastc = c;
    c = fgetc(f);
  }
  if (feof(f))
    error->all(FLERR,"fix contact/map: Could not read map file.");

  if (nx_ <= 0)
    error->all(FLERR,"fix contact/map: Something wrong: nx_ <= 0");


  /*
   * for now assume a square map
   */
  allocate_map(map_data, nx_, nx_);


  /*
   * read map data
   */
  rewind(f);
  for (int i = 0; i < nx_; i++) {
    for (int j = 0; j < ny_; j++) {
      double tmp;
      if (fscanf(f, "%lf", &tmp) != 1)
	error->all(FLERR,"fix contact/map: Error read map file.");

      int di = i - xlo_loc_;
      int dj = j - ylo_loc_;

      if (di < 0)     di += nx_;
      if (di >= nx_)  di -= nx_;
      if (dj < 0)     dj += ny_;
      if (dj >= ny_)  dj -= ny_;

      if (di >= 0 && di < nx_loc_ && dj >= 0 && dj < ny_loc_) {
	map_data[di][dj] = tmp;
      }
    }
  }

  fclose(f);

  qx_ = nx_/domain->xprd;
  qy_ = ny_/domain->yprd;
}

/* ----------------------------------------------------------------------
   write a file containing a surface map
------------------------------------------------------------------------- */

void FixWallMap::write_map(char *fn, double **map_data)
{
  if (comm->nprocs > 1)
    error->all(FLERR,"FixContacMap::write_map only works with single processor "
	       "runs.");

  FILE *f = fopen(fn, "w");

  if (!f) {
    char errstr[1024];
    sprintf(errstr, "Could not open contact map file '%s' for writing.", fn);
    error->all(FLERR,errstr);
  }

  for (int i = 0; i < nx_; i++) {
    for (int j = 0; j < ny_; j++) {
      fprintf(f, "%f ", map_data[i][j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}

/* ----------------------------------------------------------------------
   filter map for a certain low-wavelength cutoff
------------------------------------------------------------------------- */

void FixWallMap::filter_map(double **map_data)
{
#ifndef NO_LAMMPS
  // debug
  //  write_map("tmp1.matrix", map_data);

  if (comm->me == 0) {
    if (screen)
      fprintf(screen, "\nFiltering map for short wavelength cutoff %f.\n",
	      lambda_);
    if (logfile)
      fprintf(logfile, "\nFiltering map for short wavelength cutoff %f.\n",
	      lambda_);
  }

  double dqx  = 1.0/domain->xprd;
  double dqy  = 1.0/domain->yprd;
  double cutq = 1.0/lambda_;

  // allocate FFT object and temporary data buffer
  int nfft;
  FFT3d *fft = new FFT3d(lmp, world,
	       //  fast, med,                 slow
		   1,    ny_,                  nx_,
		   0, 0, ylo_loc0_, yhi_loc0_, xlo_loc0_, xhi_loc0_,
		   0, 0, ylo_loc0_, yhi_loc0_, xlo_loc0_, xhi_loc0_,
           0, 0, &nfft, 0);

  double *fft_data = (double *)
    memory->smalloc(nxy_loc0_*sizeof(FFT_DATA), "FixWallMap::fft_data");

  // fill FFT buffer (which is complex)
  int m = 0;
  for (int x = xlo_loc0_; x <= xhi_loc0_; x++) {
    for (int y = ylo_loc0_; y <= yhi_loc0_; y++) {
      fft_data[m++] = map_data[x-xlo_loc_][y-ylo_loc_];
      fft_data[m++] = 0.0;
    }
  }

  // perform the FFT!
  fft->compute(fft_data, fft_data, -1);

  // cut wavevectors
  m = 0;
  for (int x = xlo_loc0_; x <= xhi_loc0_; x++) {
    for (int y = ylo_loc0_; y <= yhi_loc0_; y++) {
      int dx = MIN(x, nx_-x);
      int dy = MIN(y, ny_-y);
      double r = sqrt((dx*dqx)*(dx*dqx) + (dy*dqy)*(dy*dqy));
      
      if (r > cutq) {
	fft_data[m  ] = 0.0;
	fft_data[m+1] = 0.0;
      }

      m += 2;
    }
  }

  // perform reverse FFT!
  fft->compute(fft_data, fft_data, 1);

  // fill map buffer
  int nxy = nx_*ny_;
  m = 0;
  for (int x = xlo_loc0_; x <= xhi_loc0_; x++) {
    for (int y = ylo_loc0_; y <= yhi_loc0_; y++) {
      map_data[x-xlo_loc_][y-ylo_loc_] = fft_data[m]/nxy;
      m += 2;
    }
  }

  // free stuff
  memory->sfree(fft_data);

  delete fft;

  // debug
  //  write_map("tmp2.matrix", map_data);
#else
  error->all(FLERR,"filter_map not available without LAMMPS FFT.");
#endif
}

/* ----------------------------------------------------------------------
   displace map such that it is aligned either with a given average or
   a given top/bottom value
------------------------------------------------------------------------- */

void FixWallMap::displace_map(double **map_data)
{
  // compute average, min, max
  double za_loc = 0.0, h0_loc = 0.0;
  double maxz_loc = map_data[0][0], minz_loc = map_data[0][0];
  for (int i = xlo_loc0_; i <= xhi_loc0_; i++) {
    for (int j = ylo_loc0_; j <= yhi_loc0_; j++) {
      double v = map_data[i-xlo_loc_][j-ylo_loc_];

      za_loc += v;
      h0_loc += v*v;

      maxz_loc = MAX(maxz_loc, v);
      minz_loc = MIN(minz_loc, v);
    }
  }

  double za, maxz, minz;

  // collect data from other procs
  MPI_Allreduce(&za_loc, &za, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&h0_loc, &h0_, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&minz_loc, &minz, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&maxz_loc, &maxz, 1, MPI_DOUBLE, MPI_MAX, world);

  h0_ = sqrt( h0_/(nx_*ny_) - (za/(nx_*ny_))*(za/(nx_*ny_)) );
  double avgh = za/(nx_*ny_);

  // remove average height and add zero height
  switch (z0_type_) {
  case BOTTOM:
    za = z0_ - minz;
    break;
  case TOP:
    za = z0_ - maxz;
    break;
  default:
    za = z0_ - za/(nx_*ny_);
    break;
  }

  for (int i = 0; i < nxy_loc_; i++) {
    map_data[0][i] += za;
  }
  avgh += za;

  /*
   * Write stuff to log file
   */
  if (comm->me == 0) {
    if (screen)
      fprintf(screen, "\nContact map has dimensions %ix%i and extends from "
	      "%f to %f with rms height %f and average %f.\n",
	      nx_, ny_, minz+za, maxz+za, h0_, avgh);
    if (logfile)
      fprintf(logfile, "\nContact map has dimensions %ix%i and extends from "
	      "%f to %f with rms height %f and average %f.\n",
	      nx_, ny_, minz+za, maxz+za, h0_, avgh);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_table2d(double x, double y, double &f,
				 double &dfdx, double &dfdy,
				 double &d2fdxdx, double &d2fdydy,
				 double &d2fdxdy)
{
  double rx = x*qx_ - xlo_loc_;
  double ry = y*qy_ - ylo_loc_;

  map_->eval(rx, ry, f, dfdx, dfdy, d2fdxdx, d2fdydy, d2fdxdy);

  dfdx *= qx_;
  dfdy *= qy_;
  d2fdxdx *= qx_*qx_;
  d2fdydy *= qy_*qy_;
  d2fdxdy *= qx_*qy_;
}

/* ---------------------------------------------------------------------- */

void FixWallMap::eval_raw_data(double x, double y, double &z,
				  double &dz_dx, double &dz_dy,
				  double &d2z_dxdx, double &d2z_dydy,
				  double &d2z_dxdy)
{
  int ix = (int) (x*qx_);
  int iy = (int) (y*qy_);
  WRAPX(ix);
  WRAPY(iy);
  z = map_data_[ix][iy];
  dz_dx = 0.0;
  dz_dy = 0.0;
  d2z_dxdx = 0.0;
  d2z_dydy = 0.0;
  d2z_dxdy = 0.0;
}
