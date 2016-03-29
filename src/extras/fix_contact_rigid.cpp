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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_contact_rigid.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixContactRigid::FixContactRigid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int iarg;
  char *endptr;
  char errstr[120];
  
  if (narg < 4) error->all(FLERR,"fix contact/rigid: Illegal command, syntax is 'fix ID group-ID contact/rigid <N> hist <np> <minp> <maxp>'");

  np = -1;
  output_forces = false;
  prefix = NULL;

  output_interval = strtol(arg[3], &endptr, 10);
  if (endptr == arg[3])
    error->all(FLERR,"fix contact/rigid: Illegal command, syntax is 'fix ID group-ID contact/rigid <N> hist <np> <minp> <maxp>'");

  // optional arguments
  iarg = 4;
  while (iarg < narg) {
    if (!strcmp(arg[iarg], "hist")) {
      if (iarg+3 >= narg)
	error->all(FLERR,"fix contact/rigid: Optional *hist* parameter requires three argument.");

      np = strtol(arg[iarg+1], &endptr, 10); 
      if (endptr == arg[iarg+1])
	error->all(FLERR,"fix contact/rigid: Number of pressure histogram points must be integer.");
      minp = strtod(arg[iarg+2], &endptr);
      if (endptr == arg[iarg+2])
	error->all(FLERR,"fix contact/rigid: Minimum pressure must be real value.");
      maxp = strtod(arg[iarg+3], &endptr);
      if (endptr == arg[iarg+3])
	error->all(FLERR,"fix contact/rigid: Maximum pressure must be real value.");

      minpsq = minp*minp;
      maxpsq = maxp*maxp;

      iarg += 3;
    }
    else if (!strcmp(arg[iarg], "output_forces")) {
      if (narg < iarg+2)
	error->all(FLERR,"fix contact/rigid: Expected prefix after output_forces keyword.");
      prefix = strdup(arg[iarg+1]);
      output_forces = true;
      iarg++;
    }
    else {
      sprintf(errstr, "fix contact/rigid: Unknown keyword '%s' encountered.",
	      arg[iarg]);
      error->all(FLERR,errstr);
    }

    iarg++;
  }

  // fix behavior

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 0;
  extvector = 1;

  // other stuff

  force_flag = 0;
  frigid_all[0] = frigid_all[1] = frigid_all[2] = 0;
  ncontact_sum = 0;
  nit = 0;
  nout = 0;

  // allocate pressure histogram

  phist_loc = NULL;
  phist_sum = NULL;
  if (np > 0) {
    memory->create(phist_loc, np+1, "fix:phist_loc");
    memory->create(phist_sum, np+1, "fix:phist_sum");
  }

  MPI_Comm_rank(world, &rank);
}

/* ---------------------------------------------------------------------- */

FixContactRigid::~FixContactRigid()
{
  if (phist_loc)
    memory->destroy(phist_loc);
  if (phist_sum)
    memory->destroy(phist_sum);
  if (prefix)
    free(prefix);
}

/* ---------------------------------------------------------------------- */

int FixContactRigid::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  nevery = 1;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactRigid::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all(FLERR,"fix contact/rigid: RESPA not yet supported.");
}

/* ---------------------------------------------------------------------- */

void FixContactRigid::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRigid::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRigid::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int tag_enable = atom->tag_enable;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  FILE *outfile;
  char fn[80];

  int k;
  double fabssq;

  /* 
   * Output forces, compute force histogram, only with frequency
   * of *output_interval*.
   */

  nit += 1;
  if (output_interval > 0 && nit >= output_interval
      && ( output_forces || phist_loc )) {

    if (phist_loc)
      memset(phist_loc, 0, (np+2)*sizeof(double));

    sprintf(fn, "%s.force.%i.%i.out", prefix, rank, nout);
    outfile = fopen(fn, "w");

    if (!atom->tag_enable)
      fprintf(outfile, "# Atom indices are local!\n");

    k = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)  k++;
    }
    fprintf(outfile, "%i\n", k);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (tag_enable) {
	  fprintf(outfile,
		  "%i %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
		  tag[i],
		  x[i][0], x[i][1], x[i][2],
		  f[i][0], f[i][1], f[i][2]);
	}
	else {
	  fprintf(outfile,
		  "%i %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
		  i,
		  x[i][0], x[i][1], x[i][2],
		  f[i][0], f[i][1], f[i][2]);
	}

	fabssq = f[i][0]*f[i][0] + f[i][0]*f[i][1] + f[i][0]*f[i][2];

	if (phist_loc) {
	  if (fabssq >= minpsq) {
	    if (fabssq < maxpsq) {
	      fabssq = sqrt(fabssq);
	      k = (fabssq-minp)*np/(maxp-minp)+1;
	      if (k < 1) k = 1;
	      if (k > np+1) k = np+1;
	      phist_loc[k]++;
	    }
	    else
	      phist_loc[np+1]++;
	  }
	  else
	    phist_loc[0]++;
	}
      }
    }

    fclose(outfile);

    if (phist_loc && phist_sum)
      MPI_Allreduce(phist_loc, phist_sum, np+2, MPI_DOUBLE, MPI_SUM, world);

    if (phist_sum && rank == 0) {
      sprintf(fn, "%s.phist.%i.out", prefix, nout);
      outfile = fopen(fn, "w");
      for (k = 0; k < np+2; k++) {
	fprintf(outfile, "%f %f\n",
		minp+(k-1+0.5)*(maxp-minp)/np, phist_sum[k]);
      }
      fclose(outfile);
    }

    nit = 0;
    nout++;
  }

  /*
   * Set forces to zero and count number of contacts.
   */

  force_flag = 0;
  ncontact_loc = 0;

  frigid[0] = frigid[1] = frigid[2] = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      frigid[0] += f[i][0];
      frigid[1] += f[i][1];
      frigid[2] += f[i][2];

      fabssq = f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
      if (fabssq > 0.0)
	ncontact_loc += 1;

      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
    }
  }

  MPI_Allreduce(&ncontact_loc, &ncontact_sum, 1, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixContactRigid::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return the contact area
------------------------------------------------------------------------- */

double FixContactRigid::compute_scalar()
{
  return ncontact_sum;
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixContactRigid::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(frigid,frigid_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return frigid_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixContactRigid::memory_usage()
{
  double bytes = 0.0;
  if (phist_loc) bytes = (np+2) * sizeof(double);
  if (phist_sum) bytes = (np+2) * sizeof(double);
  return bytes;
}
