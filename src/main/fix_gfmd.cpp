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
/* ----------------------------------------------------------------------
   Contributing authors:
     L.T. Kong, G. Bartels, C. Campana, C. Denniston and M. Muser
     
   Contact:
     Department of Applied Mathematics, University of Western Ontario
     London, ON, Canada N6A 5B7

     mmuser@uwo.ca, cdennist@uwo.ca, konglt@gmail.com

   Contributing authors:
     Lars Pastewka, Tristan A. Sharp

   Contact:
     Department of Physics and Astronomy, Johns Hopkins University,
     Baltimore, MD 21218, USA

     pas@pha.jhu.edu, tsharp@pha.jhu.edu

   Changes/extensions:
     Jan. 2011
     Transfer matrix formulation of the Green's function, in particular
     specific implementation for FCC lattices.

     Feb. 2011
     Code reorganization: separated computations of dynamical matrices
     and Green's functions from LAMMPS specifics parts.

     Mar. 2011
     Code reorganization: moved solver into a separate module;
     moved grid ids and xeq to new gfmd atom type.
     Analytical kernel for farther than nearest neighbor interactions.
     CUDA solver.

     Apr.-May 2011
     Make code more robust, in particular with respect to shift in the
     initial grid definition.

     June 2011, April 2012
     Added dynamic GFMD.

     May 2012
     Added GFMD analyzer.

     April 2013
     Integration with USER-CUDA, keep data on GPU to eliminate
     copy operations.
     Made code more memory efficient.

     Sep. 2013
     Note atom grid indexing must begin with (0,0,0) located within one 
       unit cell of the simulation box's min(x) and min(y).  
     Removed requirement in the error_checking (check_grid_indices) that 
       the atom with gid=(0,0,0) has the lowest x- and y-coord of all GF atoms.
     Dont communicate borders in fix_gfmd::setup; it's by LAMMPS in min/run.cpp
     Corrected the fix_gfmd::check_grid_indices function.
     Added in the invsurfvecec where required.
     Remove devtimer, because it required pointers to the lmp LAMMPS object.

     Mar.-Mai 2014
     Made computation of force-constants and stiffness matrix more general.
     Added support for arbitrary lattices and finite differences computation
     of force-constant matrices.

     Apr 2021
     Maintenance update for compatibility with latest LAMMPS master.
------------------------------------------------------------------------- */

#include <algorithm>
#include <cmath>
#include <iostream>

#include <fenv.h>
#include <math.h>
#include <string.h>

#include "mpi.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "fix_gfmd.h"
#include "force.h"
#include "group.h"
#include "lattice.h"
#include "linearalgebra.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include "gfmd_solver.h"
#include "surface_stiffness.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE   256

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

#define nearbyint(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define rint(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define round(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define lround(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))

#define LOC_IDX(ix, iy, iu)  ( ( (ix)*(ny_loc) + (iy) )*(nu) + (iu) )

/* ---------------------------------------------------------------------- */

FixGFMD::FixGFMD(LAMMPS *lmp, int narg, char **arg)
  : Fix(lmp, narg, arg)
{
#if 0
    /* Uncomment to enable floating-point exception */
    //    int feenableexcept();
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  int *mask = atom->mask;
  int **gid = atom->gid;
  int nlocal = atom->nlocal;

  int iarg;
  char errstr[120];

  if (!atom->gfmd_flag) {
    error->all(FLERR,"fix gfmd requires atom attributes gid and xeq; use "
               "gfmd atom style.");
  }

  if (comm->procgrid[2] > 1) {
    error->all(FLERR,"fix gfmd does not work with domain decomposition along "
               "the z-axis.");
  }

  me = comm->me;
  nprocs = comm->nprocs;

  if (narg<5) error->all(FLERR,"fix gfmd: Number of arguments < 5");
 
  iarg = 3;
  prefix = strdup(arg[iarg]);
  iarg++;

  gfmdlog = NULL;

  kernel = NULL;
  solver = NULL;

  u_xy = NULL;
  f_xy = NULL;
  f_i = NULL;
  force_is_known = NULL;

  if ((iarg+2)>narg) {
    error->all(FLERR,"fix gfmd: Insufficient command line options. "
               "Please provide keyword for stiffness expression.");
  }  

  iarg++;
  kernel = stiffness_kernel_factory(arg[iarg-1], narg, &iarg, arg,
                                    domain, force, memory, error);

  if (!kernel) {
    sprintf(errstr, "fix gfmd: Unknown analytic stiffness keyword "
            "'%s' encountered.", arg[iarg-1]);
    error->all(FLERR,errstr);
  }

  // get kernel dimension
  ndof = kernel->get_dimension();
  ndof_sq = ndof*ndof;
  //nu = kernel->get_number_of_atoms();
  nu = ndof/3;

  // to read other command line options, they are optional
  do_dump_stiffness = false;
  do_dump_greens_function = false;
  dumpr_every = 0;
  dumpq_every = 0;

  reset_xeq = false;
  reset_map = false;

  fix_center[0] = false;
  fix_center[1] = false;
  fix_center[2] = false;

  border = 0.015;

  while (iarg < narg){
    if (strcmp(arg[iarg],"fix_center") == 0) {
      if ((iarg+4)>narg)
        error->all(FLERR,"fix gfmd: fix_center: Insufficient command line "
                   "options");
      for (int i = 0; i < 3; i++) {
        iarg++;
        fix_center[i] = atoi(arg[iarg]);
      }
    }
    else if (strcmp(arg[iarg],"invariant") == 0) {
      kernel->set_invariant(true);
    }
    else if (strcmp(arg[iarg], "solver") == 0) {
      if ((iarg+2)>narg) error->all(FLERR,"fix gfmd: Insufficient command line options");
      iarg+=2;
      solver = gfmd_solver_factory(arg[iarg-1], lmp, narg, &iarg, arg);
      if (!solver) {
        sprintf(errstr, "fix gfmd: Unknown solver type '%s' "
                "encountered.", arg[iarg-1]);
        error->all(FLERR,errstr);
      }
      iarg--;
    }
    // whether/how to reset xeq based on surface lattice info from gfc
    else if (strcmp(arg[iarg],"reset_xeq") ==0){
      reset_xeq = true;
    }
    // whether/how to reset xeq based on surface lattice info from gfc
    else if (strcmp(arg[iarg],"reset_map") ==0){
      reset_map = true;
    }
    // dump stiffness matrix elements to file
    else if (!strcmp(arg[iarg], "dump_stiffness")) {
      do_dump_stiffness = true;
    }
    // dump greens functions matrix elements to file
    else if (!strcmp(arg[iarg], "dump_greens_function")) {
      do_dump_greens_function = true;
    }
    // dump interval for forces, displacements and energies (real space)
    else if (!strcmp(arg[iarg], "dumpr_every")) {
      dumpr_every = atoi(arg[++iarg]);
    }
    // dump interval for forces, displacements and energies (q-space)
    else if (!strcmp(arg[iarg], "dumpq_every")) {
      dumpq_every = atoi(arg[++iarg]);
    }
#if 0
    // specify number of iterations for Green's functions evaluation
    else if (!strcmp(arg[iarg], "gf_nlayers")) {
      gf_nlayers = atoi(arg[++iarg]);
    }
#endif
    // border multiplier for communication of atoms
    else if (!strcmp(arg[iarg], "border")) {
      iarg++;
      border = atof(arg[iarg]);
      if (border < 0.01)
        error->all(FLERR,"fix gfmd: Please provide at least 1% border.");
    }
    else {
      char str[MAXLINE];
      sprintf(str,"fix gfmd: Unknown fix arg: %s", arg[iarg]);
      error->all(FLERR,str);
    }

    iarg++;
  } // end of reading command line options

  epot = 0.0;
  epot_flag = false;
  pre_flag = false;
  grid_checked_flag = false;

  pre_time = 0.0;
  post_time = 0.0;

  // get the total number and mass (inverse) of atoms in group
  natoms   = group->count(igroup);
  masstotal  = group->mass(igroup);
  rmasstotal = 1./masstotal;
  if (natoms<1) error->all(FLERR,"fix gfmd: no atom is passed to gfmd");

  // open the log file on root
  if (me == 0){
    char str[MAXLINE];
    sprintf(str,"%s.log",prefix);
    if (logfile) {
      gfmdlog = logfile;
    }
    else {
      gfmdlog = fopen(str,"w");
    }
    if (gfmdlog == NULL){
      sprintf(str,"Cannot open log file: %s.log",prefix);
      error->one(FLERR,str);
    }
    if (fix_center[0])
      fprintf(gfmdlog, "Fixing x-center of the elastic half space.\n");
    if (fix_center[1])
      fprintf(gfmdlog, "Fixing y-center of the elastic half space.\n");
    if (fix_center[2])
      fprintf(gfmdlog, "Fixing z-center of the elastic half space.\n");
    fprintf(gfmdlog, "Using a %f %% relative thickness for the ghost shell.\n",
            border);
    fflush(gfmdlog);
  }

  if (kernel->get_height() == -1) { //TAS height == 0 means N=0 (alpha=0 only compliant layer)
    if (me == 0) {
      fprintf(gfmdlog, "Elasticity kernel returns zero height. Elastic layer "
              "will behave translationally invariantly by enforcing the q=0 sum "
              "rule.\n");
    }
    kernel->set_invariant(true);
  }

  if (gfmdlog)
    kernel->dump_info(gfmdlog);

  if (reset_xeq) {
    comp_xeq();
  }

  rnu = nu;
  if (reset_map) {
    comp_map();
    //check_grid_indices();
  }
  else {
    nx_loc = 0;
    ny_loc = 0;
    rnu = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = gid[i][0];
        int iy = gid[i][1];
        int iu = gid[i][2];

        nx_loc = std::max(nx_loc, ix);
        ny_loc = std::max(ny_loc, iy);
        rnu    = std::max(rnu,   iu);
      }
    }

    MPI_Allreduce(&nx_loc, &nx, 1, MPI_INT, MPI_MAX, world);
    MPI_Allreduce(&ny_loc, &ny, 1, MPI_INT, MPI_MAX, world);

    nx++;
    ny++;
    rnu++;

    if (rnu != nu) {
      if (screen)
        fprintf(screen, "Warning: Stiffness kernel reported %i atoms per unit "
                "cell, but there are %i atoms per unit cell in the data file "
                "for the group specified.\n", nu, rnu);
      if (gfmdlog)
        fprintf(gfmdlog, "Warning: Stiffness kernel reported %i atoms per unit "
                "cell, but there are %i atoms per unit cell in the data file "
                "for the group specified.\n", nu, rnu);
    }
  }
  nxy = nx*ny;

  if (gfmdlog) {
    fprintf(gfmdlog, "Elastic grid has dimension %i x %i with %i atoms per "
            "unit cell.\n", nx, ny, rnu);
    fflush(gfmdlog);
  }

  if (!solver) {
    // initialize default solver
    int dummy;
    solver = gfmd_solver_factory(NULL, lmp, 0, NULL, NULL);
  }
  solver->set_grid_size(nx, ny, ndof);

  xlo_loc    = solver->get_xlo_loc();
  xhi_loc    = solver->get_xhi_loc();
  ylo_loc    = solver->get_ylo_loc();
  yhi_loc    = solver->get_yhi_loc();
  nx_loc     = xhi_loc - xlo_loc + 1;
  ny_loc     = yhi_loc - ylo_loc + 1;
  nxy_loc    = mynq = mynpt = solver->get_nxy_loc();
  natoms_localsubgrid = nxy_loc*rnu;

  // output mapping info to log file
  if (gfmdlog) {
    fprintf(gfmdlog, "Elastic grid has dimension %i x %i per processor.\n",
            nx_loc, ny_loc);
    if (me == 0) {
      fprintf(gfmdlog, "%f x %f - %f x %f section of the domain is on proc "
              "0.\n", domain->sublo[0], domain->sublo[1], domain->subhi[0],
              domain->subhi[1]);
      fprintf(gfmdlog, "%i x %i - %i x %i section of the grid is on proc "
              "0.\n", xlo_loc, ylo_loc, xhi_loc, yhi_loc);
      for (int i = 1; i < nprocs; i++) {
        double buf[8];
        MPI_Status status;
        MPI_Recv(buf, 8, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
        fprintf(gfmdlog, "%f x %f - %f x %f\n section of the domain is on "
                "proc %i.\n", buf[0], buf[1], buf[2], buf[3], i);
        fprintf(gfmdlog, "%i x %i - %i x %i section of the grid is on proc "
                "%i.\n", int(buf[4]), int(buf[5]), int(buf[6]), int(buf[7]),
                i);
      }
    }
    fprintf(gfmdlog, "%d atoms found in the elastic manifold total.\n", natoms);
    fprintf(gfmdlog, "Using %s FFT solver.\n", solver->get_name());
    fflush(gfmdlog);
  }
  else {
    double buf[8];
    
    buf[0] = domain->sublo[0];
    buf[1] = domain->sublo[1];
    buf[2] = domain->subhi[0];
    buf[3] = domain->subhi[1];
    buf[4] = xlo_loc;
    buf[5] = ylo_loc;
    buf[6] = xhi_loc;
    buf[7] = yhi_loc;
    MPI_Send(buf, 8, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  /*
   * enable to return original forces before they are changed
   */
  energy_global_flag = 1;
  extvector = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  comm_forward = 3;
  comm_reverse = 4;
} // end of constructor

/* ---------------------------------------------------------------------- */

FixGFMD::~FixGFMD()
{
  if (dumpr_every > 0)
    dump();

  // delete arrays allocated by new
  free(prefix);

  // dellocate compute buffers
  memory->destroy(u_xy); u_xy = NULL;
  memory->destroy(f_xy); f_xy = NULL;
  memory->destroy(f_i); f_i = NULL;

  memory->sfree(force_is_known); force_is_known = NULL;

  // destroy FFT plan and worksapce
  delete solver;
  
  // delete stiffness kernel
  delete kernel;

  if (gfmdlog) {
    fprintf(gfmdlog, "GFMD pre_force time = %f\n", pre_time);
    fprintf(gfmdlog, "GFMD post_force time = %f\n", post_time);
    fprintf(gfmdlog, "GFMD total time = %f\n", pre_time+post_time);
    if (gfmdlog != logfile)
      fclose(gfmdlog);
  }
}

/* ---------------------------------------------------------------------- */

int FixGFMD::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGFMD::init()
{
  int *mask = atom->mask;
  int **gid = atom->gid;
  int nlocal = atom->nlocal;

  // check for another GFMD fix and warn
  int count = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"gfmd") == 0) count++;
  }
  if (count > 1 && me == 0) error->warning(FLERR,"More than one fix gfmd"); // it is allowed

  solver->init();

  solver->set_kernel(kernel);

  if (do_dump_stiffness) {
    solver->dump_stiffness();
  }
  if (do_dump_greens_function) {
    solver->dump_greens_function();
  }

  /*
   * Compute center of mass of the positions given by xeq
   */

  get_cm_from_xeq(xeqcm);

  if (gfmdlog) {
    fprintf(gfmdlog, "Center-of-mass of the elastic layer is %f, %f, %f.\n",
            xeqcm[0], xeqcm[1], xeqcm[2]);
    fflush(gfmdlog);
  }


  /*
   * If we assume 3 % max strain, we need local box size * 3%/2 spacing
   * on each side. Note: This is just an estimate, but I can't think of
   * anything better.
   */
  const mat<double> cell(3, kernel->get_cell(), true);
  double cutx = std::max(abs(cell[0][0]), abs(cell[1][0]))
    *(domain->subhi[0]-domain->sublo[1])*border;
  double cuty = std::max(abs(cell[0][1]), abs(cell[1][1]))
    *(domain->subhi[1]-domain->sublo[1])*border;
  comm->cutghostuser = std::max(comm->cutghostuser, std::max(cutx, cuty));

  if (gfmdlog) {
    fprintf(gfmdlog, "Absolute thickness of the ghost shell is %f (%f, %f).\n",
            comm->cutghostuser, cutx, cuty);
    fflush(gfmdlog);
  }


  /*
   * initialize remaining stuff
   */
  // FIXME! force_is_known is nmax in size, should be able to reduce to nall
  nmax = atom->nmax;
  xshift = 0;
  yshift = 0;

  /*
   * Check grid indices
   */
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int ix = gid[i][0];
      int iy = gid[i][1];
      int iu = gid[i][2];

      if (ix < 0 || ix >= nx) {
          char errstr[1024];
          sprintf(errstr, "Grid index x (=%i) out of range [0,%i] for atom %i.",
                  ix, nx-1, i);
          error->one(FLERR, errstr);
      }
      if (iy < 0 || iy >= ny) {
          char errstr[1024];
          sprintf(errstr, "Grid index y (=%i) out of range [0,%i] for atom %i.",
                  iy, ny-1, i);
          error->one(FLERR, errstr);
      }
      if (iu < 0 || iu >= rnu) {
          char errstr[1024];
          sprintf(errstr, "Grid index u (=%i) out of range [0,%i] for atom %i.",
                  iu, rnu-1, i);
          error->one(FLERR, errstr);
      }

    }
  }


  // allocate memory.  Number of GF atoms will never change.

  // u_xy is stored only locally.  It is double** (fix_gfmd.h) and so create() template in memory.h
  if (u_xy == NULL) memory->create(u_xy, rnu*ndof/nu, nxy_loc, "FixGFMD:u_xy");
  // f_xy has global dimension, because we need a border for reverse
  // communication (FIXME! reduce dimension to (nx+2)*(ny+2))
  if (f_xy == NULL) memory->create(f_xy, rnu*ndof/nu, nxy_loc, "FixGFMD:f_xy");
  if (f_i == NULL) memory->create(f_i, nmax, 3, "FixGFMD:f_i");
  // auxiliary buffer
  if (force_is_known == NULL) force_is_known = (bool *) memory->smalloc(nmax*sizeof(bool),
                                            "FixGFMD:force_is_known");
} //end ::init()

/* ---------------------------------------------------------------------- */

void FixGFMD::setup(int vflag)
{
  // Test whether all grid indices have been properly defined in accordance
  // with the grid defined by xlo*, yhi* etc
  // If no internal error occurs or atoms lost, grid indices should only 
  // need to be checked at init, ie when reading in the datafile from the user.
  // Setup is called each time a new run or minimize command is called.
  check_grid_indices();    

  if (!pre_flag)
    pre_force(vflag);
  post_force(vflag);

}

/* ---------------------------------------------------------------------- */

void FixGFMD::min_setup(int vflag)
{
  setup(vflag);
}

/* ---------------------------------------------------------------------- */
// pre_force updates the atom->gid's after atoms may have been exchanged 
// between procs, then fills in the displacement-on-the-grid array, 
// then updates the shift variables.
// 
// In detail, pre_force performs the following:
// Learn how many grid points the gf layer c.o.m. has moved from equilbrium (cur_xshift)
// Compare the shift to the previously-recorded shift (dxshift = xshift - cut_xshift)
// dxshift needs to be taken into account in gid, so ix=ix_from_gid(i)-dxshift
// If dx(or y)shift !=0 wrap ix to box as needed, and update gid.
// Fill in the displacement-of-grid array u_xy
// Update xshift to be the just-calculated cur_xshift
void FixGFMD::pre_force(int vflag)
{
  double **x = atom->x;
  double **xeq = atom->xeq;
  double **v = atom->v;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double xprd_half = domain->xprd_half;
  double yprd_half = domain->yprd_half;
  double *h = domain->h;

  int *mask = atom->mask;
  int *image = atom->image;
  int **gid = atom->gid;
  int *gflag = atom->gflag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int natoms_cur;

  /* --- */

  last_timestep = update->ntimestep;

  double previous_time = MPI_Wtime();

  if (nmax < atom->nmax) {
    nmax = atom->nmax;
    memory->grow(f_i, nmax, 3, "FixGFMD::f_i");
    force_is_known = (bool *) 
      memory->srealloc(force_is_known, nmax*sizeof(bool),
                       "FixGFMD:force_is_known");
  }

  double dxcm[3];
  group->xcm(igroup, masstotal, dxcm);

  dxcm[0] -= xeqcm[0];
  dxcm[1] -= xeqcm[1];
  dxcm[2] -= xeqcm[2];

  /*
   * Compute the overall shift of the GFMD layer in terms of lattice
   * coordinates. We need to shift our internal indices by this value.
   * Only tested for square grids.
   */
  const mat<double> invcell(3, kernel->get_invcell(), true);
  int cur_xshift = (int) nearbyint(invcell[0][0]*dxcm[0] +
                                   invcell[0][1]*dxcm[1] +
                                   invcell[0][2]*dxcm[2]);
  int cur_yshift = (int) nearbyint(invcell[1][0]*dxcm[0] +
                                   invcell[1][1]*dxcm[1] +
                                   invcell[1][2]*dxcm[2]);

  int dxshift = xshift - cur_xshift;
  int dyshift = yshift - cur_yshift;

  /*
  if (dxshift != 0 || dyshift != 0) {
    if (me == 0 && screen)
      fprintf(screen, "shifting: dxshift = %i, dyshift = %i\n", dxshift,
              dyshift);
  }
  */

  epot = 0.0;
  epot_flag = false;

  /*
   * get displacements U(r) on local proc
   */
  ngfmd_loc = 0;
  natoms_cur = 0;
  if (domain->triclinic == 0) {  // for orthogonal lattice
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int ix = gid[i][0] - dxshift;
        int iy = gid[i][1] - dyshift;
        int iu = gid[i][2];
        int flag = gflag[i];

        if (dxshift != 0 || dyshift != 0) {
          // wrap to box if ix or iy has shifted
          while (ix >= nx)  ix -= nx;
          while (ix < 0)    ix += nx;
          while (iy >= ny)  iy -= ny;
          while (iy < 0)    iy += ny;

          // store wrapped grid id for unpacking/communication
          gid[i][0] = ix;
          gid[i][1] = iy;
          gid[i][2] = iu;
          gflag[i] = flag;
        }

        ix -= xlo_loc;
        iy -= ylo_loc;

        double uz = x[i][2] - xeq[i][2];

        if (ix >= 0 && ix < nx_loc && iy >= 0 && iy < ny_loc) {

          // y is the fast index
          int iloc = ix*ny_loc + iy;
          int idof = 3*iu;

          // In cases where the system is 1-proc-wide in a periodic dim, 
          // I believe we can have both a local and ghost copy of an atom.
          // If we prefer i<nlocal, could use if ((i<nlocal) || !setyet(gid(i)))
          double ux = x[i][0] - xeq[i][0];
          double uy = x[i][1] - xeq[i][1];

          // wrap to box
          while (ux >  xprd_half)  ux -= xprd;
          while (ux < -xprd_half)  ux += xprd;
          while (uy >  yprd_half)  uy -= yprd;
          while (uy < -yprd_half)  uy += yprd;

          u_xy[idof  ][iloc] = ux;
          u_xy[idof+1][iloc] = uy;
          u_xy[idof+2][iloc] = uz;


#if 1
          /* 
           * Velocites of ghost atoms are undefined so we need to copy them.
           * Velocities are needed for the dynamic GFMD only.
           */
          if (i < nlocal) {
            f_xy[idof  ][iloc] = v[i][0];
            f_xy[idof+1][iloc] = v[i][1];
            f_xy[idof+2][iloc] = v[i][2];
          }
#endif

          natoms_cur++;
        }

        if (i < nlocal) {
          ngfmd_loc++;
        }
      }
    }
  }
  else {                      // for non-orthogonal lattice
    error->all(FLERR,"The non-orthogonal GF case has never been tested. "
               "Fix it!");

    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int ix = gid[i][0] - xlo_loc;
        int iy = gid[i][1] - ylo_loc;

        // wrap to box
        while (ix > nx)  ix -= nx;
        while (ix < 0)   ix += nx;
        while (iy > ny)  iy -= ny;
        while (iy < 0)   iy += ny;

        if (ix >= 0 && ix < nx_loc && iy >= 0 && iy < ny_loc) {

          int iu = gid[i][2];

          // y is the fast index
          int iloc = ix*ny_loc + iy;
          int idof = 3*iu;

          int xbox = (image[i] & 1023) - 512;
          int ybox = (image[i] >> 10 & 1023) - 512;
          int zbox = (image[i] >> 20) - 512;

          u_xy[idof  ][iloc] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox
            - xeq[i][0];
          u_xy[idof+1][iloc] = x[i][1] + h[1]*ybox + h[3]*zbox
            - xeq[i][1];
          u_xy[idof+2][iloc] = x[i][2] + h[2]*zbox
            - xeq[i][2];

          natoms_cur++;
        }
      }
    }
  }  // end of orthogonal/ non statements 

  if (natoms_cur != natoms_localsubgrid) {
    char errstr[1024];
    sprintf(errstr, "natoms_cur != natoms_localsubgrid (%i != %i)",
            natoms_cur, natoms_localsubgrid);
  }

  xshift = cur_xshift;
  yshift = cur_yshift;

  solver->pre_force(u_xy, f_xy); // virtual function in gfmd_solver.h

  pre_time += MPI_Wtime() - previous_time;

  pre_flag = true;
}


/* ---------------------------------------------------------------------- */
// post_force call solver->post_force to get the [restoring] forces
// restoring force put into list format and added to LAMMPS' force
void FixGFMD::post_force(int vflag)
{
  double **f = atom->f;

  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  int nall   = nlocal + atom->nghost;

  double previous_time = MPI_Wtime();

  if (!pre_flag)
    error->one(FLERR,"FixGFMD: pre_force needs to be called before "
               "post_force");

  if (dumpq_every > 0 && update->ntimestep % dumpq_every == 0) {
    char dump_prefix[1024];
    sprintf(dump_prefix, "%s.%li", prefix, update->ntimestep);
    epot += solver->post_force(u_xy, f_xy, dump_prefix);
  }
  else 
    epot += solver->post_force(u_xy, f_xy, NULL);

  /*
   * Copy f_xy from grid back to force list f_i
   */
  grid_to_list();
  fsum_flag = false;

  /*
   * Add forces from list to atoms
   */

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += f_i[i][0];
      f[i][1] += f_i[i][1];
      f[i][2] += f_i[i][2];  
    }
  }

  /*
   * Dump displacements, forces and energies
   */
  if (dumpr_every > 0 && update->ntimestep % dumpr_every == 0)
    dump();

  /*
   * Remove force on the center of mass
   */
  if (fix_center[0] || fix_center[1] || fix_center[2]) {
    MPI_Allreduce(fsum_loc, fsum, 3, MPI_DOUBLE, MPI_SUM,
                  world);
    fsum_flag = true;

    for (int idim = 0; idim < 3; idim++) {
      if (fix_center[idim]) {
        double df = fsum[idim]/natoms;

        for (int i = 0; i < nall; i++) {
          if (mask[i] & groupbit) {
            f[i][idim] -= df;
          }
        }
      }
    }
  }

  /*
   * Gather total energy on all processors
   */
  double tmp;
  MPI_Allreduce(&epot, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  epot = tmp;

  epot_flag = true;
  pre_flag  = false;

  post_time += MPI_Wtime() - previous_time;

} // end post_force(int)

/* ----------------------------------------------------------------------
 * copy force data from the grid f_xy back to the f_i list that stores
 * the forces per atom. f_i also contains forces on ghost atoms while
 * f_xy exclusively stores information on grid atoms. As implemented, all
 * lammps-local atoms and all gfmd-local atoms will have a force known.
 * --------------------------------------------------------------------*/

void FixGFMD::grid_to_list()
{
  int **gid = atom->gid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  /*
   * add elastic force to local atoms
   */
  fsum_loc[0] = fsum_loc[1] = fsum_loc[2] = 0.0;
  for (int i = 0; i < nall; i++) {
    force_is_known[i] = false;
  }

  for (int i = 0; i < nall; i++) { 
    if (mask[i] & groupbit) {
      int ix = gid[i][0];
      int iy = gid[i][1];

      // has this atom been computed locally?
      if (ix >= xlo_loc && ix <= xhi_loc && iy >= ylo_loc && iy <= yhi_loc) {

          ix -= xlo_loc;
          iy -= ylo_loc;

        int iu = gid[i][2];

        // y is the fast index
        int iloc = ix*ny_loc + iy;
        int idof = 3*iu;

        // GFMD forces
        double fx = f_xy[idof  ][iloc];
        double fy = f_xy[idof+1][iloc];
        double fz = f_xy[idof+2][iloc];
        
        // we store the GFMD force separately
        f_i[i][0] = fx;
        f_i[i][1] = fy;
        f_i[i][2] = fz; 
        // f_i is thus set for all "gfmd_local" atoms, whether lmp_local or lmp_ghost atoms.
        // f_i will be set in reverse_comm_fix for atoms that are gfmd_nonlocal and lmp_local.
        // A force_is_known error-check below ensures each local gfmd atom is assigned a gfmd-force.

        // sum all forces on GFMD layer
        if (i < nlocal) {
          fsum_loc[0] += fx;
          fsum_loc[1] += fy;
          fsum_loc[2] += fz;
        }

        if (force_is_known[i]) {
          error->one(FLERR,"Duplicate find.");
        }
        force_is_known[i] = true;
      }
    }
  }

  /*
   * communicate other procs' ghost data back to this processor
   */
  // f_i is here set for all atoms that are both "gfmd_nonlocal" and ghosts of other procs.
  comm->reverse_comm_fix(this);

  // f_i is now set for all lmp_local atoms
  // as long as all local atoms were either gfmd_local 
  // or else reverse-communicated and gfmd_nonlocal.
  // This all assumes gfmd_local can be found in the list of local+ghost atoms. 
  // And that lmp_local atoms can be found in either gfmd_local atoms or reverse-comm'd border region.

  // If forces on all of my proc's ghost atoms are additionally required
  // comm->forward_comm_fix(this); is also available, though may need to test it first.

  /*
   * Check that all local atoms have a gfmd force assigned
   */

  double **x   = atom->x;
  double **xeq = atom->xeq;
  int *tag     = atom->tag;

  double xcm[3];
  group->xcm(igroup, masstotal, xcm);

  double *u0 = solver->get_u0();

  char loststr[1024];
  strcpy(loststr, "");
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (!force_is_known[i]) {
        char newstr[1024];

        int ix = gid[i][0];
        int iy = gid[i][1];
        int iu = gid[i][2];

        sprintf(newstr, "%i (%i %i %i - %f %f %f),\n ", tag[i],
                ix, iy, iu, x[i][0]-xeq[i][0], x[i][1]-xeq[i][1],
                x[i][2]-xeq[i][2]);
        if (strlen(loststr) + strlen(newstr) < 1024-2)
          strcat(loststr, newstr);
      }
    }
  }

  if (strlen(loststr) >= 2) {
    char errstr[2048];
    loststr[strlen(loststr)-2] = 0;

    sprintf(errstr, "fix gfmd: Not all local atoms received a GFMD force.\n Are "
            "the initial grid positions set correctly in the input file?\n "
            "If system distored, consider increasing comm border distance.\n "
            "Overall GFMD layer center of mass is %f %f %f. q=0\n "
            "displacement is %f %f %f. The atoms that were lost are \n %s.",
            xcm[0], xcm[1], xcm[2], u0[0], u0[1],
            u0[2], loststr);
    error->one(FLERR,errstr);
  }

}

/* ---------------------------------------------------------------------- */

void FixGFMD::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGFMD::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixGFMD::pack_forward_comm(int n, int *list, double *buf, int pbc_flag,
                               int *pbc)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];

    buf[m++] = f_i[j][0];
    buf[m++] = f_i[j][1];
    buf[m++] = f_i[j][2];
  }

  return 3;
}

/* ---------------------------------------------------------------------- */

void FixGFMD::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    // This would require selecting for only those atoms outside this 
    // proc's subgrid range.  if ((ix < xlo_loc) || (ix > xhi_loc))...
    f_i[i][0] = buf[m++];
    f_i[i][1] = buf[m++];
    f_i[i][2] = buf[m++];
    force_is_known[i] = true;
  }
}

/* ---------------------------------------------------------------------- */

int FixGFMD::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    // If force_is_known, ie was a local_gfmd atom here, then 
    // communicate the GFMD force back
    if (force_is_known[i]) {
      buf[m++] = i-first;
      buf[m++] = f_i[i][0];
      buf[m++] = f_i[i][1];
      buf[m++] = f_i[i][2];
    }
  }
  // termination
  buf[m++] = n;
  return m;
}

/* ---------------------------------------------------------------------- */
// Fill in f_i and fsum_loc for the GFMD-non-local atoms that are our ghost atoms.
// Developer's note: LAMMPS does not normally share the forces on LAMMPS-local atoms.
// buf has [int i, fx tot on atom j, fy, fz,   int i+1, ...]
// where int i is such that list[i]== atom_j_index 
// "j < nlocal" means "owned by this processor"
// LAMMPS will give us a list of atom local-indices n long, and their forces in buf
void FixGFMD::unpack_reverse_comm(int n, int *list, double *buf)
{
  int **gid = atom->gid;
  int nlocal = atom->nlocal;

  int m = 0;

  int i = (int) buf[m++];
  if (n > 0) {
    while (i < n) {
      int j = list[i];
      int ix = gid[j][0];
      int iy = gid[j][1];

      // Only accept this gfmd force for atom i if it wasn't already calculated here.
      // pack_reverse_comm only sends atoms that were gfmd_local on the neighboring proc.
      // However, if data is sent to self, then may receive atom that was calculated here.
      if (ix < xlo_loc || ix > xhi_loc || iy < ylo_loc || iy > yhi_loc) {
        int iu = gid[j][2];

        // #ifdef GFMD_DEBUG
        if (force_is_known[j]) {
          if (abs(f_i[j][0] - buf[m]) > 1e-12 ||
              abs(f_i[j][1] - buf[m+1]) > 1e-12 ||
              abs(f_i[j][2] - buf[m+2]) > 1e-12)
            error->one(FLERR,"Duplicate find with different forces.");
        }
        //#endif

        force_is_known[j] = true;

        f_i[j][0] = buf[m++];
        f_i[j][1] = buf[m++];
        f_i[j][2] = buf[m++];

        if (j < nlocal) {
          fsum_loc[0] += f_i[j][0];
          fsum_loc[1] += f_i[j][1];
          fsum_loc[2] += f_i[j][2];
        }
      }
      else {
        m += 3;
      }
      
      // This is the other proc's local-index for this atom, minus "first". 
      // Remind me, why compare i against n (instead of m)? 
      // n is the number of ghost atoms on this proc.
      i = (int) buf[m++]; 
    }
  }
}

/* ----------------------------------------------------------------------
   effective potential energy of Green's functions solid
------------------------------------------------------------------------- */

double FixGFMD::compute_scalar()
{
  if (!epot_flag) {
    if (!pre_flag)
      pre_force(0);
    post_force(0);
  }
  return epot;
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixGFMD::compute_vector(int n)
{
  if (!fsum_flag){
    MPI_Allreduce(fsum_loc, fsum, 3, MPI_DOUBLE, MPI_SUM, world);
    fsum_flag = true;
  }
  return fsum[n];
}

/* ----------------------------------------------------------------------
 * private method, to extract the equilibrium positions of atoms in GFMD
 * layer from the initial configuration.
 * For restart run, the initial configuration frequently does not
 * corresponds to the equilibrium position, one should provide the 
 * informaiton by file, instead.
 * --------------------------------------------------------------------*/

void FixGFMD::comp_xeq()
{
  double **x   = atom->x;
  double **xeq = atom->xeq;
  int *mask    = atom->mask;
  int *image   = atom->image;
  int nlocal   = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], xeq[i]);
    }
  }

  if (gfmdlog) {
    fprintf(gfmdlog,"Original positions of GF atoms determined from initial "
            "configuration!\n");
  }
}

/* ----------------------------------------------------------------------
 * private method, to calculate the mapping info for surface atoms.
 * --------------------------------------------------------------------*/

void FixGFMD::comp_map()
{
  double **xeq = atom->xeq;
  int *type = atom->type;
  int *mask = atom->mask;
  int **gid = atom->gid;
  int *gflag = atom->gflag;
  int nlocal = atom->nlocal;

  int my_nu;

  /* --- */
  
  // reset to zero
  for (int i = 0; i < nlocal; ++i) {
    gid[i][0] = 0;
    gid[i][1] = 0;
    gid[i][2] = 0;
    gflag[i] = 0;
  }

  // invcell is currently in 1/nn-units but should be in 1/distance units!
  // Working assumption: invcell is (100 010 001) (surface_stiffness.cpp)
  const mat<double,3> invcell(kernel->get_invcell());
  
  // get FFT dimensions
  nx = int(invcell[0][0]*domain->xprd+invcell[0][1]*domain->yprd+0.1);
  ny = int(invcell[1][0]*domain->xprd+invcell[1][1]*domain->yprd+0.1);
  if (nx <= 0 || ny <= 0) {
      char errstr[1024];
      sprintf(errstr, "Something is wrong. My guess for the dimension of the "
              "GFMD grid is %i x %i.", nx, ny);
      error->all(FLERR,errstr);
  }
  my_nu = natoms / (nx*ny);

  if (gfmdlog) {
    fprintf(gfmdlog,"I think this system has a lattice of %i x %i "
            "with %i atoms per unit cell.\n", nx, ny, my_nu);
  }

  // sanity checks
  if (nx < 1 || nx > natoms || ny < 1 || ny > natoms) {
    char errstr[1024];
    sprintf(errstr, "FFT dimensions (nx = %i, ny = %i) are likely wrong.",
            nx, ny);
    error->all(FLERR,errstr);
  }

  if (my_nu != nu) {
    char errstr[1024];
    sprintf(errstr, "Number of atoms per unit cell from stiffness kernel "
            "(= %i) and input file (= %i) differ.", nu, my_nu);
    error->all(FLERR,errstr);
  }

  if (nx*ny*nu != natoms) {
    char errstr[1024];
    sprintf(errstr, "Number of atoms (= %i) differs from nx*ny*nu (= %i).",
            natoms, nx*ny*nu);
    error->all(FLERR,errstr);
  }

  // find topmost GF atom
  double xeq_of_highest_gf_layer = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xeq_of_highest_gf_layer = std::max(xeq_of_highest_gf_layer, xeq[i][2]);
    }
  }

  // calculate the x,y-mapping info
  int j = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int ix, iy;
 
      double scaled_xeq[3];  // equilibirum position in grid units

      // to get the fractional coordination of atom i with the basis of
      // surface vectors
      double tmp[3];
      tmp[0] = xeq[i][0] - domain->boxlo[0];
      tmp[1] = xeq[i][1] - domain->boxlo[1];
      tmp[2] = xeq[i][2] - xeq_of_highest_gf_layer;
      mat_mul_vec(3, invcell.const_data(), tmp, scaled_xeq);

      ix = (int)floor(scaled_xeq[0]+1e-6);
      iy = (int)floor(scaled_xeq[1]+1e-6);

      ix %= nx;
      iy %= ny;
      ix = (ix < 0)? (ix+nx):ix;
      iy = (iy < 0)? (iy+ny):iy;

      gid[i][0] = ix;
      gid[i][1] = iy;
      gid[i][2] = 0;
      gflag[i] = 0;

      j++;
    }
  }

  xlo_loc = nx*(domain->sublo[0]-domain->boxlo[0])/
      (domain->boxhi[0]-domain->boxlo[0]);
  xhi_loc = nx*(domain->subhi[0]-domain->boxlo[0])/
      (domain->boxhi[0]-domain->boxlo[0])-1;
  ylo_loc = ny*(domain->sublo[1]-domain->boxlo[1])/
      (domain->boxhi[1]-domain->boxlo[1]);
  yhi_loc = ny*(domain->subhi[1]-domain->boxlo[1])/
      (domain->boxhi[1]-domain->boxlo[1])-1;

  // Sanity check.
  if (xlo_loc == xhi_loc) {
    error->one(FLERR,"Lowest and highest x-index of the grid are indentical. "
               "This is likely a bug. Contact the developer.");
  }
  if (ylo_loc == yhi_loc) {
    error->one(FLERR,"Lowest and highest y-index of the grid are indentical. "
               "This is likely a bug. Contact the developer.");
  }

  nx_loc = xhi_loc-xlo_loc+1;
  ny_loc = yhi_loc-ylo_loc+1;
  
  if (nu > 1) {
    // We have a multi-atom unit cell. Calculate the u-mapping info by
    // comparing to the unit cell as obtained from the stiffness kernel.

    // Collect a single unit cell from input file.
    // We use the one in the middle of the simulation domain because that one
    // is for sure on this MPI process.
    int refx = (xlo_loc+xhi_loc)/2;
    int refy = (ylo_loc+yhi_loc)/2;
    int k = 0;
    // *ref_positions* contains a 3x3 section of unit cells. This is needed
    // because atoms can be wrongly assigned to the neighboring cell in the
    // procedure that determines ix and iy above.
    double ref_positions[3*3*3*nu];
    int ref_types[3*3*nu];

    // We first identify the central cell.
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = gid[i][0];
        int iy = gid[i][1];

        if (ix == refx && iy == refy) {
          // Sanity check.
          if (k >= nu) {
            error->one(FLERR,"Unit cell index *k* exceeded number of unit cell "
                       "atoms when identifying the lattice. This is likely a "
                       "bug. Contact the developer.");
          }
          ref_positions[3*k] = xeq[i][0];
          ref_positions[3*k+1] = xeq[i][1];
          ref_positions[3*k+2] = xeq[i][2];
          ref_types[k] = type[i];
          k++;
        }
      }
    }

    // Next, we identify the 8 neighboring cells.
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = gid[i][0];
        int iy = gid[i][1];

        if (std::abs(ix-refx) <= 1 && std::abs(iy-refy) <= 1 &&
            !(ix == refx && iy == refy)) {
          // Sanity check.
          if (k >= 3*3*nu) {
            error->one(FLERR,"Unit cell index *k* exceeded number of unit cell "
                       "atoms when identifying the lattice. This is likely a "
                       "bug. Contact the developer.");
          }
          ref_positions[3*k] = xeq[i][0];
          ref_positions[3*k+1] = xeq[i][1];
          ref_positions[3*k+2] = xeq[i][2];
          ref_types[k] = type[i];
          k++;
        }
      }
    }

    // Sanity check.
    if (k != 3*3*nu) {
      char errstr[1024];
      sprintf(errstr, "Could not find all unit cell atoms. "
              "(k = %i, 3*3*nu = %i)", k, 3*3*nu);
      error->one(FLERR,errstr);
    }

    // Mapping from positions stored in *ref_positions* to the atom indices for
    // the respective unit cell.
    int indices[3*3*nu];

    if (!kernel->get_mapping(3*3*nu, ref_positions, indices, ref_types)) {
        error->all(FLERR,"Could not find mapping from current input file to "
                   "unit cell reported by the stiffness kernel. Are the "
                   "lattice type and lattice constant correct?");
    }

    // Copy the mapping to a basis.
    double basis[3*nu];
    int basis_types[nu];
    for (int i = 0; i < 3*3*nu; i++) {
        int k = indices[i];
        if (k >= 0) {
            basis[3*k] = ref_positions[3*i];
            basis[3*k+1] = ref_positions[3*i+1];
            basis[3*k+2] = ref_positions[3*i+2];
            basis_types[k] = ref_types[i];
        }
    }

    // Dump basis information to file.
    if (gfmdlog) {
        fprintf(gfmdlog,"Found mapping from input file to unit cell reported "
                "by the stiffness kernel. Prototype:\n");
        for (int i = 0; i < nu; i++) {
            fprintf(gfmdlog,"  u=%i --- ( %f, %f, %f ) --- atom type %i\n", i,
                    basis[3*i], basis[3*i+1], basis[3*i+2], basis_types[i]);
        }
    }

    // Get information about manybody kernels.
    bool manybody = kernel->is_manybody();
    int flags[nu];
    std::fill(flags, flags+nu, -1);

    // Now we use the mapping to assign grid-indices.
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = -1, iy = -1, iu = -1;

        for (int j = 0; j < nu; j++) {
          double d[3] = { xeq[i][0]-basis[3*j],
                          xeq[i][1]-basis[3*j+1],
                          xeq[i][2]-basis[3*j+2] };
          vec<double,3> xd = invcell.dot(d);
          // Check that this atom is an integer number of unit cells away from
          // the reference.
          if (std::abs(xd[0]-nearbyint(xd[0])) < 1e-3 && 
              std::abs(xd[1]-nearbyint(xd[1])) < 1e-3 && 
              std::abs(xd[2]) < 1e-3) {
            ix = refx+nearbyint(xd[0]);
            iy = refy+nearbyint(xd[1]);
            iu = j;
            if (type[i] != basis_types[j]) {
              char errstr[1024];
              sprintf(errstr,"Atom with grid index %i x %i x %i has type %i, "
                      "but should have type %i.", ix, iy, iu, type[i],
                      basis_types[j]);
              error->one(FLERR,errstr);
            }
          }
        }

        if (iu < 0 || iu >= nu) {
          char errstr[1024];
          sprintf(errstr, "Could not find u-index for local atom %i with "
                  "equilibrium position %e %e %e.", i, xeq[i][0], xeq[i][1],
                  xeq[i][2]);
          error->one(FLERR,errstr);
        }

        ix %= nx;
        iy %= ny;
        ix = (ix < 0)? (ix+nx):ix;
        iy = (iy < 0)? (iy+ny):iy;

        // Tell the interatomic potential to compute energy for just half of
        // the lattice planes. (flag=1 means do not compute energy for that
        // lattice plane.)
        int flag = 1;
        if (manybody) {
            flag = iu >= nu/2 ? 1 : 0;
        }
        flags[iu] = flag;
        gid[i][0] = ix;
        gid[i][1] = iy;
        gid[i][2] = iu;
        gflag[i] = flag;
      }
    }

    if (gfmdlog) {
      fprintf(gfmdlog,"Flags per layer:");
      for (int iu = 0; iu < nu; iu++) fprintf(gfmdlog," %i",flags[iu]);
      fprintf(gfmdlog,"\n");
    }
  }

#if 0
  for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
          int ix = gid[i][0];
          int iy = gid[i][1];
          int iu = gid[i][2];
          int flag = gflag[i];

          std::cout << vec<double>(3, xeq[i]) << ", " << ix << ", " << iy << ", " << iu << ", " << flag << std::endl;
      }
  }
#endif

  if (gfmdlog) {
    fprintf(gfmdlog,"Mapping info computed from initial configuration.\n");
  }
}

/* ----------------------------------------------------------------------
 * Test whether grid indices are missing or doubly defined
 *
 * Each atom has a grid (ix,iy,iu) equilibrium (x0,y0,z0) and
 * actual current location (x,y,z). 
 * 
 * Use several atoms' gid's to find a (0,0) ref location for each gfmd
 * plane (iu).  Then imagine a checkerboard laid down, and find if atoms 
 * equilibrium positions lie within the cell given by their gid.
 *
 * The atoms' xeq don't have to lie within their correct squares; a global 
 * shift of the gids may have occurred, and check that the shift index
 * is the same for all atoms, shift_of_ref_atom is global.
 *
 * Also keep track of which gids have been seen in the local atoms in 
 * the "exists" array.  If a gid is missing, that means this proc is 
 * missing an atom with the correct xeq.
 *
 * Although the grid squares may be shifted relative to the equilibrium
 * grid, we check here that it is not distorted out of the squares at all.
 * An error is thrown if each atom in the grid x_i was not paired with 
 * one in atom->x.
 * --------------------------------------------------------------------*/

void FixGFMD::check_grid_indices()
{
  double **xeq = atom->xeq;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int **gid = atom->gid;
  int *gflag = atom->gflag;
  int nall = atom->nlocal + atom->nghost;

  double xprd = domain->xprd;
  double yprd = domain->yprd;

  int nx2 = nx/2;
  int ny2 = ny/2;
  int xshift_of_ref_atom;
  int yshift_of_ref_atom;

#ifdef HAVE_C99
  double xreflo[rnu];
  double yreflo[rnu];
#else
  double *xreflo = new double[rnu];
  double *yreflo = new double[rnu];
#endif

  /* --- */

  const mat<double> invcell(3, kernel->get_invcell(), true);

  natoms_localsubgrid = nx_loc*ny_loc*rnu;
  int *exists = (int *) memory->smalloc(natoms_localsubgrid*sizeof(int), 
                                        "FixGFMD::comp_map::exists");

  /*
   * Determine the reference x,yreflo locations for each gfmd plane iu.
   */
  if (me == 0) {

    /*
     * look for first atom in layer rnu-1 (typically minimum-z layer) on proc 0
     */
    int i = 0;
    while (!(mask[i] & groupbit && gid[i][2] == rnu-1) &&
           i < nall)  i++;
    if (i >= nall) {
      char errstr[256];
      sprintf(errstr, "fix gfmd: No atoms for GFMD layer %d found.", rnu-1);
      error->one(FLERR,errstr);
    }

    /*
     * If we lay checkerboard down from the domain corner, in what
     * grid square is the equilibrium position of this random atom i
     */
    // To avoid truncation error when comparing atom locations to grid cells.
    // If atom is within 1e-3 of cell's upper-right-hand corner, belongs to the next grid cell.
    const double epsilon = 1e-3 * xprd / nx; 
    int ix0eq = int(floor((xeq[i][0] + epsilon - domain->boxlo[0]) * nx / xprd)); 
    int iy0eq = int(floor((xeq[i][1] + epsilon - domain->boxlo[1]) * ny / yprd));

    /*
     * For a reference point, extrapolate to where we would expect 
     * to find the identical atom in the lowest leftmost grid square, 
     * that is, in the (0,0) grid square.  This will be the reflo variables.
     */
    xreflo[rnu-1] = xeq[i][0] - ix0eq * xprd/nx;
    yreflo[rnu-1] = xeq[i][1] - iy0eq * yprd/ny;

    /*
     * Get (post any shift) grid information of atom i
     */
    int ix0 = gid[i][0];
    int iy0 = gid[i][1];

    if (gfmdlog) {
      fprintf(gfmdlog, "Reference zero for layer %d is %f, %f.\n", rnu-1,
              xreflo[rnu-1], yreflo[rnu-1]);
    }

    for (int iu = 0; iu < rnu-1; iu++) {
      /*
       * look for the atom in layer iu that has the gid=ix0, iy0
       */
      int i = 0;
      while (!(mask[i] & groupbit &&
               gid[i][0] == ix0 &&
               gid[i][1] == iy0 &&
               gid[i][2] == iu) &&
             i < nall)  i++;
      if (i >= nall) {
        char errstr[1024];
        sprintf(errstr, "fix gfmd: Atom with grid index %i x %i x %i "
                "not found.", ix0, iy0, iu);
        error->one(FLERR,errstr);
      }

      /*
       * extrapolate that ix0,iy0 back (ix0eq squares) to get this layer's reference point
       */
      xreflo[iu] = xeq[i][0] - ix0eq * xprd/nx;
      yreflo[iu] = xeq[i][1] - iy0eq * yprd/ny;

      // If the grid index for this atom doesnt match this atom's equilbrium, the reflo will 
      // be off for this layer and we will catch the error below
      if (gfmdlog) {
        fprintf(gfmdlog, "Reference zero for layer %i is %f, %f.\n",
                iu, xreflo[iu], yreflo[iu]);
      }
    }

    // We will see if the gid matches the equilibrium positioms, up to a global shift.

    // Nominally the atom->gid's would only be shifted if the shift variable were non-zero.
    // However, each time a run or minimization is called, the shift variables are re-zeroed,
    // And this check_grid_indices function is called at setup.

    // Thus relearn the shift applied to these gids based on our reference atom and 
    // its equilibrium position.  Then we will set shift to this.
    xshift_of_ref_atom = ix0 - ix0eq;
    yshift_of_ref_atom = iy0 - iy0eq;

    while (xshift_of_ref_atom < 0)   xshift_of_ref_atom += nx;
    while (xshift_of_ref_atom >= nx) xshift_of_ref_atom -= nx;
    while (yshift_of_ref_atom < 0)   yshift_of_ref_atom += ny;
    while (yshift_of_ref_atom >= ny) yshift_of_ref_atom -= ny;

  }

  MPI_Bcast(xreflo, rnu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(yreflo, rnu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xshift_of_ref_atom, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&yshift_of_ref_atom, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /*
   * We will update the gids and the shift variables
   * Compute the overall shift of the GFMD layer in terms of lattice
   * coordinates. We need to shift our gid internal indices by this value.
   */
  double dxcm[3];
  group->xcm(igroup, masstotal, dxcm);

  dxcm[0] -= xeqcm[0]; // c.o.m. of actual loc <minus> c.o.m. of equilib loc
  dxcm[1] -= xeqcm[1];
  dxcm[2] -= xeqcm[2];
  int cur_xshift = (int) nearbyint(invcell[0][0]*dxcm[0] +
                                   invcell[0][1]*dxcm[1] +
                                   invcell[0][2]*dxcm[2]);
  int cur_yshift = (int) nearbyint(invcell[1][0]*dxcm[0] +
                                   invcell[1][1]*dxcm[1] +
                                   invcell[1][2]*dxcm[2]);
  /*
   * xshift and cur_xshift are the same quantity, it's just xshift may be stale.
   * For example, if there was a shift "right before the restart file was written."
   */
  int dxshift = xshift - cur_xshift;
  int dyshift = yshift - cur_yshift;


  /*
   * the ghost shell is 2 lattice constants
   */

  /*
   * mark all grid positions as nonexistent
   */
  for (int i = 0; i < natoms_localsubgrid; i++)
    exists[i] = -1;

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      int ix = gid[i][0]; 
      int iy = gid[i][1];
      int iu = gid[i][2];
      int flag = gflag[i];

      /*
       * Compute difference between stored grid index (ix) and what we'd expect from
       * equilibrium positions (ixeq=xeq-xreflo). This should be identical for all 
       * atoms (and equal to xshift_of_ref_atom).
       */
      int ixeq = lround((xeq[i][0] - xreflo[iu]) * nx / xprd);
      int iyeq = lround((xeq[i][1] - yreflo[iu]) * ny / yprd);
      int dix = ix - ixeq;
      int diy = iy - iyeq;
      
      while (dix < 0)    dix += nx;
      while (dix >= nx)  dix -= nx;
      while (diy < 0)    diy += ny;
      while (diy >= ny)  diy -= ny;

      // This does not check that the xshift variable is correctly set.
      // Cannot compare against cur_xshift since gid hasnt been updated yet
      // if atoms were just displaced.
      // Cannot compare against xshift since xshift is reinitialized with the fix 
      // but gids atom_vec properties) can be shifted by previous runs/minimizations.
      if (dix != xshift_of_ref_atom || diy != yshift_of_ref_atom) {
        char errstr[1024];

        sprintf(errstr, "fix gfmd: Atom %i has grid index %i x %i x %i. "
                "And equilibrium position %f %f %f which "
                "falls into grid index %i x %i. Grid shift is hence "
                "%i x %i. But shift in COM from COM_eq says it should be %i x %i.",
                tag[i], ix, iy, iu, xeq[i][0], xeq[i][1], xeq[i][2], 
                ixeq, iyeq, dix, diy, xshift_of_ref_atom, yshift_of_ref_atom);
        error->one(FLERR,errstr);
      }

      // Update the xshift variable and gids, with the cur_xshift
      // First update the gids
      if (xshift_of_ref_atom != cur_xshift || yshift_of_ref_atom != cur_yshift) {

        ix -= (xshift_of_ref_atom - cur_xshift);
        iy -= (yshift_of_ref_atom - cur_yshift);
        // wrap to box
        while (ix >= nx)  ix -= nx;
        while (ix < 0)    ix += nx;
        while (iy >= ny)  iy -= ny;
        while (iy < 0)    iy += ny;
      
        // store wrapped grid id
        gid[i][0] = ix;
        gid[i][1] = iy;
        gid[i][2] = iu;
        gflag[i] = flag;
        // xshift updated below, out of the loop
      }

      //printf("me = %i, ix = %i, iy = %i, nlocal = %i, nghost = %i\n", me, ix, iy, atom->nlocal, atom->nghost);

      // Make sure we can fill in the whole exists[nx,ny,nz] grid with atom tags
      // from the proc's list of local and ghost atoms, without conflicting tags.
      if (ix >= xlo_loc && ix <= xhi_loc &&
          iy >= ylo_loc && iy <= yhi_loc) {

        int j = LOC_IDX(ix-xlo_loc, iy-ylo_loc, iu);

        if (exists[j] != -1 && exists[j] != tag[i]) {
          char errstr[1024];
          sprintf(errstr, "fix gfmd: Atom %i has same grid index %i x %i x %i "
                  "as atom %i, proc %i (gfmd grid %i->%i %i->%i) reports.",
                  tag[i], ix, iy, iu, exists[j], me, xlo_loc, xhi_loc, ylo_loc, yhi_loc);
          error->one(FLERR,errstr);
        } else exists[j] = tag[i];      
      }
    }
  }

  // check if an atom is missing from the exists grid
  char errstr[10240];
  bool atom_missing=0;
  strcpy(errstr, "Atom(s) missing:\n"); 
  for (int ix = 0; ix < nx_loc; ix++) {
    for (int iy = 0; iy < ny_loc; iy++) {
      for (int iu = 0; iu < rnu; iu++) {
        int j  = LOC_IDX(ix, iy, iu);
        char newstr[1024];
        if (exists[j] < 0) {
          atom_missing=1;
          sprintf(newstr, " index = %i, grid x,y,z = "
                  "%i %i %i (MPI proc. %i) \n", j, xlo_loc+ix, ylo_loc+iy, iu, me);
          if (strlen(errstr) + strlen(newstr) < 10240-2)
            strcat(errstr, newstr);
        }
      }
    }
  }
  if (atom_missing) error->one(FLERR,errstr);

  memory->sfree(exists);

  xshift = cur_xshift;
  yshift = cur_yshift;

  grid_checked_flag = true;

#ifndef HAVE_C99
  delete [] xreflo;
  delete [] yreflo;
#endif
}

/* ----------------------------------------------------------------------
 * private method to compute the center-of-mass coordinates of GFMD atoms
 * at equilibrium positions
 * --------------------------------------------------------------------*/
 
void FixGFMD::get_cm_from_xeq(double *xcm)
{
  int *mask     = atom->mask;
  int *type     = atom->type;
  double *mass  = atom->mass;
  double *rmass = atom->rmass;
  double **xeq  = atom->xeq;
  int nlocal    = atom->nlocal;

  double xcm_loc[3];

  xcm_loc[0] = 0.0;
  xcm_loc[1] = 0.0;
  xcm_loc[2] = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        xcm_loc[0] += xeq[i][0] * mass[type[i]];
        xcm_loc[1] += xeq[i][1] * mass[type[i]];
        xcm_loc[2] += xeq[i][2] * mass[type[i]];
      }
    }
  }
  else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        xcm_loc[0] += xeq[i][0] * rmass[i];
        xcm_loc[1] += xeq[i][1] * rmass[i];
        xcm_loc[2] += xeq[i][2] * rmass[i];
      }
    }
  }

  MPI_Allreduce(xcm_loc, xcm, 3, MPI_DOUBLE, MPI_SUM, world);

  xcm[0] *= rmasstotal;
  xcm[1] *= rmasstotal;
  xcm[2] *= rmasstotal;
}

/* ---------------------------------------------------------------------- */

double FixGFMD::memory_usage()
{
  double bytes = solver->memory_usage();

  // u_xy
  bytes += nxy_loc*ndof*sizeof(double);
  // f_xy
  bytes += nxy_loc*ndof*sizeof(double);
  // f_i
  bytes += nmax*3*sizeof(double);
  // force_is_known
  bytes += nmax*sizeof(bool);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixGFMD::dump()
{
  if (nprocs > 1)
    error->all(FLERR,"fix gfmd: Can only dump from single processor run.");

  char dump_prefix[1024], fn[1024];
  sprintf(dump_prefix, "%s.%li", prefix, last_timestep);

#ifdef HAVE_C99
  FILE *fu[ndof], *ff[ndof];
#else
  FILE *fu[MAX_NDOF], *ff[MAX_NDOF];
#endif

  for (int idof = 0; idof < ndof; idof++) {
    sprintf(fn, "%s.r.u%i.out", dump_prefix, idof);
    fu[idof] = fopen(fn, "w");
    sprintf(fn, "%s.r.f%i.out", dump_prefix, idof);
    ff[idof] = fopen(fn, "w");
  }
  sprintf(fn, "%s.r.uP.out", dump_prefix);
  FILE *fuP = fopen(fn, "w");
  sprintf(fn, "%s.r.fP.out", dump_prefix);
  FILE *ffP = fopen(fn, "w");
  sprintf(fn, "%s.r.e.out", dump_prefix);
  FILE *fe = fopen(fn, "w");

  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      int iloc = ix*ny_loc + iy;

      double uP = 0.0, fP = 0.0, e = 0.0;
      for (int idof = 0; idof < ndof; idof++) {
        fprintf(fu[idof], " %20.10e ", u_xy[idof][iloc]);
        fprintf(ff[idof], " %20.10e ", f_xy[idof][iloc]);

        uP += u_xy[idof][iloc]*u_xy[idof][iloc];
        fP += f_xy[idof][iloc]*f_xy[idof][iloc];

        e += u_xy[idof][iloc]*f_xy[idof][iloc];
      }
      fprintf(fuP, " %20.10e ", uP);
      fprintf(ffP, " %20.10e ", fP);
      fprintf(fe, " %20.10e ", e);
    }
    for (int idof = 0; idof < ndof; idof++) {
      fputc('\n', fu[idof]);
      fputc('\n', ff[idof]);
    }
    fputc('\n', fuP);
    fputc('\n', ffP);
    fputc('\n', fe);
  }

  for (int idof = 0; idof < ndof; idof++) {
    fclose(fu[idof]);
    fclose(ff[idof]);
  }
  fclose(fuP);
  fclose(ffP);
  fclose(fe);
}

