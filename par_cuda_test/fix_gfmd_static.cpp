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

     Lars Pastewka
     Department of Physics and Astronomy, Johns Hopkins University,
     Baltimore, MD 21218, USA
     pas@pha.jhu.edu
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>

#include "mpi.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "fix_gfmd_static.h"
#include "force.h"
#include "group.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "linearalgebra.h"

#include "gfmd_solver.h"
#include "gfmd_stiffness.h"

using namespace LAMMPS_NS;


#define MAXLINE   256

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/*
 * Assuming 32 bit integers:
 * This limits us to GFMD sizes of 16384 x 16384 and kernel sizes of 16
 * (5 atoms per unit cell) which should be okay for any conceivable
 * situation on a 32 bit machine
 */
#define POW2_IDX(ix, iy, iu) ( (ix << 18) + (iy << 4) + iu )
#define CONT_IDX(ix, iy, iu)  ( ( (ix)*(ny) + (iy) )*(nu) + (iu) )
#define LOC_IDX(ix, iy, iu)  ( ( (ix)*(ny_loc) + (iy) )*(nu) + (iu) )

#define IX_FROM_POW2_IDX(idx)  ( idx >> 18 )
#define IY_FROM_POW2_IDX(idx)  ( ( idx >> 4 ) & 16383 )
#define IU_FROM_POW2_IDX(idx)  ( idx & 15 )

#define IX_FROM_CONT_IDX(idx)  ( (idx) / ( (ny) * (nu) ) )
#define IY_FROM_CONT_IDX(idx)  ( ( (idx) / (nu) ) % (ny) )
#define IU_FROM_CONT_IDX(idx)  ( (idx) % (nu) )

#define POW2_TO_CONT_IDX(idx)  CONT_IDX(IX_FROM_POW2_IDX(idx), IY_FROM_POW2_IDX(idx), IU_FROM_POW2_IDX(idx))

/* ---------------------------------------------------------------------- */

FixGFMDStatic::FixGFMDStatic(LAMMPS *lmp, int narg, char **arg)
  : Fix(lmp, narg, arg)
{
  int *mask  = atom->mask;
  int *gid   = atom->gid;
  int nlocal = atom->nlocal;

  int iarg;
  char errstr[120];
  char *phi_fn;
  char *solver_class_str;
  StiffnessKernel *phi_kernel;

  if (!atom->gfmd_flag) {
    error->all("fix gfmd/static requires atom attributes gid and xeq; use "
	       "gfmd atom style");
  }

  me = comm->me;
  nprocs = comm->nprocs;
  //  MPI_Comm_rank(world,&me);
  //  MPI_Comm_size(world,&nprocs);
 
  if (narg<5) error->all("fix gfmd/static: Number of arguments < 5");
 
  iarg = 3;
  prefix = strdup(arg[iarg]);
  iarg++;

  gfmdlog = NULL;

  phi_fn = NULL;
  phi_kernel = NULL;

  if ((iarg+1)>narg) {
    error->all("fix gfmd/static: Insufficient command line options. "
	       "Please provide keyword for stiffness expression.");
  }

  if (!strcmp(arg[iarg], "file")) {
    /* Read stiffness coefficients from file */
    if ((iarg+2)>narg) error->all("fix gfmd/static: Filename expected.");
    phi_fn = strdup(arg[iarg]);
    iarg += 2;
  }
  else {
    iarg++;
    phi_kernel = stiffness_kernel_factory(arg[iarg-1], narg, &iarg, arg,
					  domain, memory, error);

    if (!phi_kernel) {
      sprintf(errstr, "fix gfmd/static: Unknown analytic stiffness keyword "
	      "'%s' encountered.", arg[iarg-1]);
      error->all(errstr);
    }
  }

  // to read other command line options, they are optional
  noutfor = 0;
  do_dump_stiffness = false;
  solver_class_str = NULL;

  reset_xeq = false;
  reset_map = false;

  have_su = false;
  have_sv = false;

  fix_center[0] = false;
  fix_center[1] = false;
  fix_center[2] = false;

  while (iarg < narg){
    // surface vector U. if not given, will be determined from lattice info
    if (strcmp(arg[iarg],"su") == 0){
      if (iarg+3 > narg) error->all("fix gfmd/static: Insufficient command line options.");
      surfvec[0][0] = atof(arg[++iarg]);
      surfvec[0][1] = atof(arg[++iarg]);
      have_su = true;
    }
    // surfactor vector V. if not given for 3D, will be determined from lattice info
    else if (strcmp(arg[iarg],"sv") == 0){
      if (iarg+3 > narg)
	error->all("fix gfmd/static: Insufficient command line options.");
      surfvec[1][0] = atof(arg[++iarg]);
      surfvec[1][1] = atof(arg[++iarg]);
      have_sv = true;
    }
    else if (strcmp(arg[iarg],"fix_center") == 0) {
      if ((iarg+4)>narg)
	error->all("fix gfmd/static: fix_center: Insufficient command line "
		   "options");
      for (int i = 0; i < 3; i++) {
	iarg++;
	fix_center[i] = atoi(arg[iarg]);
      }
    }
    else if (strcmp(arg[iarg], "fft") == 0) {
      if ((iarg+2)>narg) error->all("fix gfmd/static: Insufficient command line options");
      iarg++;
      solver_class_str = strdup(arg[iarg]);
    }
    // frequency to output elastic force
    else if (strcmp(arg[iarg],"output") ==0){
      if (iarg+2 > narg) error->all("fix gfmd/static: Insufficient command line options.");
      noutfor = atoi(arg[++iarg]);
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
    // allow nonzero phi0 (i.e. finite stiffness)
    else {
      char str[MAXLINE];
      sprintf(str,"fix gfmd/static: Unknown fix arg: %s", arg[iarg]);
      error->all(str);
    }

    iarg++;
  } // end of reading command line options

  epot = 0.0;
  epot_flag = false;
  pre_flag = false;

  pre_time = 0.0;
  post_time = 0.0;

  sysdim  = domain->dimension; // find the system and surface dimension
  surfdim = sysdim-1;
  if (sysdim == 2){
    surfvec[1][0] = 0.;
    surfvec[1][1] = 1.;
    have_su = true;
  }
  ndof = phi_kernel->get_dimension();
  ndof_sq = ndof*ndof;
  nu = ndof/sysdim;

  // get the total number and mass (inverse) of atoms in group
  natoms   = group->count(igroup);
  masstotal  = group->mass(igroup);
  rmasstotal = 1./masstotal;
  if (natoms<1) error->all("fix gfmd/static: no atom is passed to gfmd");

  // open the log file on root
  if (me == 0){
    char str[MAXLINE];
    sprintf(str,"%s.log",prefix);
    gfmdlog = fopen(str,"w");
    if (gfmdlog == NULL){
      sprintf(str,"Cannot open log file: %s.log",prefix);
      error->one(str);
    }
    fprintf(gfmdlog, "GFMD initialization...\n");
    fprintf(gfmdlog, "Number of processors used    : %d\n", nprocs);
    if (phi_kernel) {
      fprintf(gfmdlog, "Analytic GF expression       : %s\n",
	      phi_kernel->getName());
    }
    if (fix_center[0])
      fprintf(gfmdlog, "Fixing x-center of GFMDlayer\n");
    if (fix_center[1])
      fprintf(gfmdlog, "Fixing y-center of GFMDlayer\n");
    if (fix_center[2])
      fprintf(gfmdlog, "Fixing z-center of GFMDlayer\n");
    fflush(gfmdlog);
  }

  if (reset_xeq) {
    comp_xeq();
  }

  if (reset_map) {
    comp_map(nu, &nx, &ny);
  }
  else {
    nx_loc = 0;
    ny_loc = 0;
    nu     = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	int j = gid[i];

	int ix = IX_FROM_POW2_IDX(j);
	int iy = IY_FROM_POW2_IDX(j);
	int iu = IU_FROM_POW2_IDX(j);

	nx_loc = max(nx_loc, ix);
	ny_loc = max(ny_loc, iy);
	nu     = max(nu,     iu);
      }
    }

    MPI_Allreduce(&nx_loc, &nx, 1, MPI_INT, MPI_MAX, world);
    MPI_Allreduce(&ny_loc, &ny, 1, MPI_INT, MPI_MAX, world);

    nx++;
    ny++;
    nu++;
  }

  nxy = nx*ny;

  if (gfmdlog) {
    fprintf(gfmdlog, "GF grid has dimension        : %i x %i\n", nx, ny);
    fprintf(gfmdlog, "# of atoms per unit cell     : %d\n", nu);
    fflush(gfmdlog);
  }

  // reset surface vectors based on box info
  surfvec[0][0] = domain->h[0]/double(nx);
  surfvec[0][1] = domain->h[5]/double(nx);
  surfvec[1][0] = 0.;
  surfvec[1][1] = domain->h[1]/double(ny);
  if (sysdim == 2) surfvec[1][1] = 1.;

  solver = gfmd_solver_factory(solver_class_str, lmp, nx, ny, ndof);
  if (!solver) {
    sprintf(errstr, "fix gfmd/static: Unknown FFT type '%s' encountered.",
	    solver_class_str);
    error->all(errstr);
  }
  if (solver_class_str)
    free(solver_class_str);

  xlo_loc    = solver->get_xlo_loc();
  xhi_loc    = solver->get_xhi_loc();
  ylo_loc    = solver->get_ylo_loc();
  yhi_loc    = solver->get_yhi_loc();
  nx_loc     = xhi_loc - xlo_loc + 1;
  ny_loc     = yhi_loc - ylo_loc + 1;
  nxy_loc    = mynq = mynpt = solver->get_nxy_loc();
  natoms_loc = nxy_loc*nu;

  // u_xy is stored only locally
  u_xy  = memory->create_2d_double_array(ndof, nxy_loc, "fix_gfmd:u_xy");
  // f_xy has global dimension, because we need a border for reverse
  // communication (FIXME! reduce dimension to (nx+2)*(ny+2))
  f_xy  = memory->create_2d_double_array(ndof, nxy,     "fix_gfmd:f_xy");

  // output mapping info to log file
  if (gfmdlog) {
    fprintf(gfmdlog, "Grid dimension per processor : %i x %i\n",
	    nx_loc, ny_loc);
    fprintf(gfmdlog, "Total # of atoms in manifold : %d\n", natoms);
    fprintf(gfmdlog, "Surface vector along U       : [ %lg, %lg ]\n",
	    surfvec[0][0],surfvec[0][1]);
    fprintf(gfmdlog, "Surface vector along V       : [ %lg, %lg ]\n",
	    surfvec[1][0],surfvec[1][1]);
    fprintf(gfmdlog, "Using FFT solver             : %s\n", solver->get_name());
    fflush(gfmdlog);
  }

  /*
   * allocate memory and get phi_q, only store those relevant to local proc
   */
  phi_q = solver->create_complex_operator_buffer("FixGFMDStatic:phi_q");
  // allocating remaining working arrays;
  // MAX(1,.. is used to avoid MPI buffer error

  if (!phi_kernel) {
    // read phi_q from binary file from FixGFC run
    readphi(phi_fn);
    free(phi_fn);
  }
  else {
    // Calculate analytic phi_q for simple cubic lattice
    phi_analytic(phi_kernel);
    delete phi_kernel;
  }

  if (do_dump_stiffness) {
    dump_stiffness();
  }

  // divide phi_q by (nx*ny) to reduce float operation after FFT
  double inv_nxny = 1./double(nx*ny);
  for (int idq=0; idq<mynq; idq++) {
    for (int idim=0; idim<ndof_sq; idim++) {
      // premptively divide for later conven.
      phi_q[idq][idim] *= inv_nxny;
    }
  }

  solver->set_operator(phi_q);


  /*
   * Compute center of mass of the positions given by xeq
   */

  get_cm_from_xeq(xeqcm);

  if (gfmdlog) {
    fprintf(gfmdlog, "Center-of-mass of GF layer   : %f %f %f\n",
	    xeqcm[0], xeqcm[1], xeqcm[2]);
    fflush(gfmdlog);
  }


  /*
   * We need at least a lattice constant thick layer of ghost atoms
   */
  double cutx = MAX(abs(surfvec[0][0]), abs(surfvec[1][0]));
  double cuty = MAX(abs(surfvec[0][1]), abs(surfvec[1][1]));
  comm->cutghostuser = MAX(comm->cutghostuser, MAX(cutx, cuty));

  if (gfmdlog) {
    fprintf(gfmdlog, "Thickness of ghost shell     : %f\n",
	    comm->cutghostuser);
    fflush(gfmdlog);
  }


  /*
   * initialize remaining stuff
   */
  nmax = atom->nmax;
  force_is_known = (bool *) memory->smalloc(nmax*sizeof(bool),
					    "FixGFMDStatic:force_is_known");
  xshift = 0;
  yshift = 0;


  /*
   * enable to return original forces before they are changed
   */
  extvector = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = sysdim;
  global_freq = 1;
  extscalar = 1;
} // end of constructor

/* ---------------------------------------------------------------------- */

FixGFMDStatic::~FixGFMDStatic()
{
  // always write the elastic force at the final step
  if (noutfor<=0 || update->ntimestep%nevery!=0) end_of_step();

  // delete arrays allocated by new
  free(prefix);
  memory->sfree(force_is_known);

  // delete arrays grow by memory->..
  memory->destroy_2d_double_array(u_xy);
  memory->destroy_2d_double_array(f_xy);
  solver->destroy_complex_operator_buffer(phi_q);

  // destroy FFT plan and worksapce
  delete solver;

  if (gfmdlog) {
    fprintf(gfmdlog, "GFMD pre_force time   = %f\n", pre_time);
    fprintf(gfmdlog, "GFMD post_force time  = %f\n", post_time);
    fprintf(gfmdlog, "GFMD total time       = %f\n", pre_time+post_time);
    fclose(gfmdlog);
  }
}

/* ---------------------------------------------------------------------- */

int FixGFMDStatic::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  if (noutfor>0){
    mask  |= END_OF_STEP;
    nevery = noutfor;
  }
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::init()
{
  int count = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"gfmd") == 0) count++;
  }
  if (count > 1 && me == 0) error->warning("More than one fix gfmd/static"); // it is allowed

} //end ::init()

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::setup(int vflag)
{
  if (!pre_flag)
    pre_force(vflag);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::min_setup(int vflag)
{
  if (!pre_flag)
    pre_force(vflag);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::pre_force(int vflag)
{
  double **x   = atom->x;
  double **xeq = atom->xeq;

  double xprd      = domain->xprd;
  double yprd      = domain->yprd;
  double zprd      = domain->zprd;
  double xprd_half = domain->xprd_half;
  double yprd_half = domain->yprd_half;
  double zprd_half = domain->zprd_half;
  double *h        = domain->h;

  int *mask  = atom->mask;
  int *image = atom->image;
  int *gid   = atom->gid;
  int nlocal = atom->nlocal;
  int nall   = nlocal + atom->nghost;

  double dxcm[3];

  int natoms_cur;
  int cur_xshift, cur_yshift;
  int dxshift, dyshift;

  /* --- */

  double previous_time = MPI_Wtime();

  if (nmax < atom->nmax) {
    nmax = atom->nmax;
    force_is_known = (bool *) 
      memory->srealloc(force_is_known, nmax*sizeof(bool),
		       "FixGFMDStatic:force_is_known");
  }

  group->xcm(igroup, masstotal, dxcm);

  dxcm[0] -= xeqcm[0];
  dxcm[1] -= xeqcm[1];
  dxcm[2] -= xeqcm[2];


  /*
   * Compute the overall shift of the GFMD layer in terms of lattice
   * coordinates. We need to shift our internal indices by this value.
   * FIXME! this should be the inverse surface vector
   */
  cur_xshift = (int) nearbyint(surfvec[0][0]*dxcm[0] + surfvec[0][1]*dxcm[1]);
  cur_yshift = (int) nearbyint(surfvec[1][0]*dxcm[0] + surfvec[1][1]*dxcm[1]);

  dxshift = xshift - cur_xshift;
  dyshift = yshift - cur_yshift;

  /*
   * get U(r) on local proc
   */
  ngfmd_loc = 0;
  natoms_cur = 0;
  if (domain->triclinic == 0) {  // for orthogonal lattice
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int j = gid[i];

	int ix = IX_FROM_POW2_IDX(j) - dxshift;
	int iy = IY_FROM_POW2_IDX(j) - dyshift;
	int iu = IU_FROM_POW2_IDX(j);

	if (dxshift != 0 || dyshift != 0) {
	  // wrap to box
	  while (ix >= nx)  ix -= nx;
	  while (ix < 0)    ix += nx;
	  while (iy >= ny)  iy -= ny;
	  while (iy < 0)    iy += ny;

	  // store wrapped grid id for unpacking/communication
	  gid[i] = POW2_IDX(ix, iy, iu);
	}

	ix -= xlo_loc;
	iy -= ylo_loc;

	if (ix >= 0 && ix < nx_loc && iy >= 0 && iy < ny_loc) {

	  // y is the fast index
	  int iloc = ix*ny_loc + iy;
	  int idof = 3*iu;

	  double ux = x[i][0] - xeq[i][0] - dxcm[0];
	  double uy = x[i][1] - xeq[i][1] - dxcm[1];
	  double uz = x[i][2] - xeq[i][2] - dxcm[2];

	  // wrap to box
	  while (ux >  xprd_half)  ux -= xprd;
	  while (ux < -xprd_half)  ux += xprd;
	  while (uy >  yprd_half)  uy -= yprd;
	  while (uy < -yprd_half)  uy += yprd;
	  while (uz >  zprd_half)  uz -= zprd;
	  while (uz < -zprd_half)  uz += zprd;

	  u_xy[idof  ][iloc] = ux;
	  u_xy[idof+1][iloc] = uy;
	  u_xy[idof+2][iloc] = uz;

	  natoms_cur++;
	}

	if (i < nlocal) {
	  ngfmd_loc++;
	}
      }
    }
  }
  else {                      // for non-orthogonal lattice
    error->all("The non-orthogonal GF case has never been tested. Fix it!");

    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int j = gid[i];

	int ix = IX_FROM_POW2_IDX(j) - xlo_loc;
	int iy = IY_FROM_POW2_IDX(j) - ylo_loc;

	// wrap to box
	while (ix > nx)  ix -= nx;
	while (ix < 0)   ix += nx;
	while (iy > ny)  iy -= ny;
	while (iy < 0)   iy += ny;

	if (ix >= 0 && ix < nx_loc && iy >= 0 && iy < ny_loc) {

	  int iu = IU_FROM_POW2_IDX(j);

	  // y is the fast index
	  int iloc = ix*ny_loc + iy;
	  int idof = 3*iu;

	  int xbox = (image[i] & 1023) - 512;
	  int ybox = (image[i] >> 10 & 1023) - 512;
	  int zbox = (image[i] >> 20) - 512;

	  u_xy[idof  ][iloc] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox
	    - dxcm[0] - xeq[i][0];
	  u_xy[idof+1][iloc] = x[i][1] + h[1]*ybox + h[3]*zbox
	    - dxcm[1] - xeq[i][1];
	  u_xy[idof+2][iloc] = x[i][2] + h[2]*zbox
	    - dxcm[2] - xeq[i][2];

	  natoms_cur++;
	}
      }
    }
  }  // TAS end of orthogonal/ non statements 

  if (natoms_cur != natoms_loc) {
    char errstr[1024];
    sprintf(errstr, "fix gfmd/static: natoms_cur != natoms_loc "
	    "(%i != %i)", natoms_cur, natoms_loc);
  }

  xshift = cur_xshift;
  yshift = cur_yshift;

  solver->detach(u_xy, f_xy);

  pre_time += MPI_Wtime() - previous_time;

  pre_flag = true;
}


/* ---------------------------------------------------------------------- */

void FixGFMDStatic::post_force(int vflag)
{
  double **f   = atom->f;

  int *mask  = atom->mask;
  int *gid   = atom->gid;
  int nlocal = atom->nlocal;
  int nall   = nlocal + atom->nghost;

  double time;

  int natoms_cur;

  double previous_time = MPI_Wtime();

  if (!pre_flag)
    error->one("FixGFMDStatic: pre_force needs to be called before "
	       "post_force");

  epot = solver->join(u_xy, f_xy);

  /*
   * add elastic force to local atoms
   */
  fsum_loc[0] = fsum_loc[1] = fsum_loc[2] = 0.0;
  fsum_flag = false;

  for (int i = 0; i < nall; i++) {
    force_is_known[i] = false;
  }

  natoms_cur = 0;

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      int j = gid[i];

      int ix = IX_FROM_POW2_IDX(j);
      int iy = IY_FROM_POW2_IDX(j);

      // has this atom been computed locally?
      if (ix >= xlo_loc && ix <= xhi_loc && iy >= ylo_loc && iy <= yhi_loc) {

	int iu = IU_FROM_POW2_IDX(j);

	// y is the fast index
	int iloc = ix*ny + iy;
	int idof = 3*iu;

	f[i][0]     += f_xy[idof  ][iloc];
	f[i][1]     += f_xy[idof+1][iloc];
	f[i][2]     += f_xy[idof+2][iloc];

	if (i < nlocal) {
	  fsum_loc[0] += f[i][0];
	  fsum_loc[1] += f[i][1];
	  fsum_loc[2] += f[i][2];
	  natoms_cur++;
	}

#ifdef GFMD_DEBUG
	if (force_is_known[i]) {
	  error->one("Duplicate find.");
	}
#endif
	force_is_known[i] = true;
      }
    }
  }


  /*
   * communicate ghost data back to this processor
   */
  natoms_loc_recv = 0;
  comm->reverse_comm_fix(this);
  natoms_cur += natoms_loc_recv;


  /*
   * Check if anything went wrong
   */
  if (natoms_cur != ngfmd_loc) {
    char errstr[1024];
    sprintf(errstr, "fix gfmd/static: natoms_cur != ngfmd_loc (%i != %i). Are "
	    "the initial grid positions set correctly in the input file?",
	    natoms_cur, ngfmd_loc);
    error->one(errstr);
  }


  /*
   * Remove force on the center of mass
   */
  if (fix_center[0] || fix_center[1] || fix_center[2]) {
    MPI_Allreduce(fsum_loc, fsum, sysdim, MPI_DOUBLE, MPI_SUM,
		  world);
    fsum_flag = true;

    for (int idim = 0; idim < sysdim; idim++) {
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

  epot_flag = true;
  pre_flag  = false;

  post_time += MPI_Wtime() - previous_time;
} // end post_force(int)

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixGFMDStatic::pack_reverse_comm(int n, int first, double *buf)
{
  int *gid = atom->gid;

  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {

    // communicate atoms for which the force is known, either because if was
    // computed locally or because it has been received earlier
    // FIXME! only send if within a single lattice constant from the border
    if (force_is_known[i]) {
      int j = gid[i];

      int ix = IX_FROM_POW2_IDX(j);
      int iy = IY_FROM_POW2_IDX(j);
      int iu = IU_FROM_POW2_IDX(j);

      // y is the fast index
      int iloc = ix*ny + iy;
      int idof = 3*iu;

      buf[m++] = i-first;
      buf[m++] = f_xy[idof  ][iloc];
      buf[m++] = f_xy[idof+1][iloc];
      buf[m++] = f_xy[idof+2][iloc];
    }

  }
  // termination
  buf[m++] = n;
  return m;
}

/* ---------------------------------------------------------------------- */

void FixGFMDStatic::unpack_reverse_comm(int n, int *list, double *buf)
{
  int *gid = atom->gid;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  int m = 0;

  int i = (int) buf[m++];
  while (i < n) {
    int j = list[i];
    int k = gid[j];
    int ix = IX_FROM_POW2_IDX(k);
    int iy = IY_FROM_POW2_IDX(k);

    // only add if this hasn't been already computed here, might happen
    // if data is sent to self
    // FIXME! only send if within a single lattice constant from the border
    if (ix < xlo_loc || ix > xhi_loc || iy < ylo_loc || iy > yhi_loc) {
      int iu = IU_FROM_POW2_IDX(k);
      int iloc = ix*ny + iy;
      int idof = 3*iu;

      f[j][0] += f_xy[idof  ][iloc] = buf[m++];
      f[j][1] += f_xy[idof+1][iloc] = buf[m++];
      f[j][2] += f_xy[idof+2][iloc] = buf[m++];

#ifdef GFMD_DEBUG
      if (force_is_known[j]) {
	error->one("Duplicate find.");
      }
#endif
      force_is_known[j] = true;

      if (j < nlocal) {
	fsum_loc[0] += f[j][0];
	fsum_loc[1] += f[j][1];
	fsum_loc[2] += f[j][2];
	natoms_loc_recv++;
      }
    }
    else {
      m += 3;
    }

    i = (int) buf[m++];
  }
}

/* ----------------------------------------------------------------------
   effective potential energy of Green's functions solid
------------------------------------------------------------------------- */

double FixGFMDStatic::compute_scalar()
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

double FixGFMDStatic::compute_vector(int n)
{
  if (!fsum_flag){
    MPI_Allreduce(fsum_loc, fsum, sysdim, MPI_DOUBLE, MPI_SUM, 
		  world);
    fsum_flag = true;
  }
  return fsum[n];
}

/* ----------------------------------------------------------------------
 * private method, to get the analytic elastic stiffness coefficients
 * for simple cubic systems
 * Loop over all wavevectors qx, qy and call the stiffness kernel to
 * get the appropriate stiffness matrix.
 * --------------------------------------------------------------------*/

void FixGFMDStatic::phi_analytic(StiffnessKernel *kernel)
{
  if (sysdim == 2){    // (1+1) D system
    error->all("FixGFMDStatic::phi_analytic: Support for 2D systems is "
	       "currently broken.");

    /*
    idq = 0;
    for (int i=xlo_loc; i<=xhi_loc; i++){
      qx = (i <= (nx/2)) ? (2.0*M_PI*i/nx) : (2.0*M_PI*(i-nx)/nx);

      K_SC_2D(qx, G_q, domain, error);
      int ndim = 0;
      for (int idim=0; idim<sysdim; idim++) {
        for (int jdim=0; jdim<sysdim; jdim++) {
	  phi_q[idq][ndim] = spring_k*G_q[ndim];
	  ndim++;
	}
      }
      idq++;
    }
    */
  }
  else {               // (2+1) D system
    int idq = 0;
    for (int i = xlo_loc; i <= xhi_loc; i++) {
      double qx = (i <= int((nx)/2)) ?
	(2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
      for (int j = ylo_loc; j <= yhi_loc; j++) {
        double qy = (j <= int((ny)/2)) ? 
	  (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

	double_complex loc_phi[ndof_sq];

	// this is where our GFunctions^-1 are imported from
        kernel->get_stiffness_matrix(qx, qy, loc_phi);
        int ndim = 0;
        for (int idim = 0; idim < ndof; idim++) {
          for (int jdim = 0; jdim < ndof; jdim++) {
	    phi_q[idq][ndim] = loc_phi[ndim];
	    ndim++;
	  }
        }
        idq++;
      }
    }
  } // end of if (sysdim ==

  //memset(phi_q[0], 0, 9*sizeof(double_complex));
}

/* ----------------------------------------------------------------------
 * private method, dump stiffness coefficients for plotting in gnuplot
 * --------------------------------------------------------------------*/

void FixGFMDStatic::dump_stiffness()
{
  int i, j, n, idim, jdim, ndim;
  FILE *freal, *fimag;
  char fn[80];

  if (me == 0) {
    ndim = 0;
    for (idim = 0; idim < ndof; idim++) {
      for (jdim = 0; jdim < ndof; jdim++) {
	sprintf(fn, "phi%i%i_real.out", idim, jdim);
	freal = fopen(fn, "w");
	sprintf(fn, "phi%i%i_imag.out", idim, jdim);
	fimag = fopen(fn, "w");

	n = 0;
	for (j = 0; j < ny; j++) {
	  for (i = xlo_loc; i <= xhi_loc; i++) {
	    fprintf(freal, "%i %i %e\n", i, j, creal(phi_q[n][ndim]));
	    fprintf(fimag, "%i %i %e\n", i, j, cimag(phi_q[n][ndim]));
	    n++;
	  }
	  fprintf(freal, "\n");
	  fprintf(fimag, "\n");
	}
	
	fclose(freal);
	fclose(fimag);

	ndim++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * private method, to extract the equilibrium positions of atoms in GFMD
 * layer from the initial configuration.
 * For restart run, the initial configuration frequently does not
 * corresponds to the equilibrium position, one should provide the 
 * informaiton by file, instead.
 * --------------------------------------------------------------------*/

void FixGFMDStatic::comp_xeq()
{
  double **x   = atom->x;
  double **xeq = atom->xeq;
  int *mask    = atom->mask;
  int *image   = atom->image;
  int nlocal   = atom->nlocal;

  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], xeq[i]);
    }
  }

  if (gfmdlog) {
    fprintf(gfmdlog,"\nOriginal positions of GF atoms determined from initial "
	    "configuration!\n");
  }
}

/* ----------------------------------------------------------------------
 * private method, to calculate the mapping info fo surface atoms.
 * For restart run, usually it is not possible to calculate this from
 * the initial configuration; instead, one would better provide the
 * mapping information by file.
 * --------------------------------------------------------------------*/

void FixGFMDStatic::comp_map(int nu, int *out_nx, int *out_ny)
{
  double **xeq = atom->xeq;
  int *mask    = atom->mask;
  int *gid     = atom->gid;
  int nlocal   = atom->nlocal;

  int my_nu, nx, ny;

  /* --- */

  if (!have_su || !have_sv) {
    // Check if lattice information is available from input file!
    if (domain->lattice == NULL){
      // get surface vectors from lattice info; orthogonal lattice is assumed
      error->all("fix gfmd/static: No lattice defined while keyword su and/or "
		 "sv not set");
    }
    if (!have_su) {
      surfvec[0][0] = domain->lattice->xlattice;
      surfvec[0][1] = 0.0;
    }
    if (!have_sv) {
      surfvec[1][0] = 0.0;
      surfvec[1][1] = domain->lattice->ylattice;
    }
  }
  // check the validity of the surface vectors read from command line
  if (fabs(surfvec[0][1]) > 0.0 && fabs(surfvec[1][0]) > 0.0)
    error->all("fix gfmd/static: Either U or V must be on the box side");
  if (surfvec[0][0] <= 0.0)
    error->all("fix gfmd/static: Surface vector U must be along the +x direction");
  if (surfvec[1][1] <= 0.0)
    error->all("fix gfmd/static: Surface vector V must point to the +y "
	       "direction");

  double invSurfV[2][2];
  for (int i=0; i<2; i++) {
    invSurfV[i][0] = surfvec[i][0];
    invSurfV[i][1] = surfvec[i][1];
  }

  // get the inverse transpose of surfvec
  GaussJordan(2,invSurfV[0], error);
  double dswap = invSurfV[0][1];
  invSurfV[0][1] = invSurfV[1][0];
  invSurfV[1][0] = dswap;

  // get FFT dimensions
  nx = int(domain->xprd*invSurfV[0][0]+0.1);
  ny = (sysdim == 2)?1:int(domain->yprd*invSurfV[1][1]+0.1);
  if (nx<1 || nx>natoms || ny<1 || ny>natoms)
    error->all("fix gfmd/static: error encountered while getting FFT "
	       "dimensions");

  my_nu = natoms / (nx*ny);
  if (my_nu != nu) {
    char errstr[1024];
    sprintf(errstr, "fix gfmd/static: number of atoms per unit cell from "
	    "stiffness kernel (= %i) and input file (= %i) differ.",
	    nu, my_nu);
    error->all(errstr);
  }

  if (nu > 1) {
    error->warning("More that one atom per unit cell: Trying to guess cell "
		   "index, but it is strongly suggested you provide the "
		   "cell indices manually.");
  }

  // now to calculate the mapping info
  int j = 0;
  xlo_loc = nx;
  xhi_loc = 0;
  ylo_loc = ny;
  yhi_loc = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int ix, iy, iu;
 
      double vx[2], vi[2];

      // relative coordinates on the surface of atom i
      vx[0] = xeq[i][0];
      vx[1] = xeq[i][1];

      // to get the fractional coordination of atom i with the basis of
      // surface vectors
      MatMulVec(2,invSurfV[0],vx,vi);

      ix = (int)floor(vi[0]+0.1); ix %= nx; ix = (ix < 0)? (ix+nx):ix;
      iy = (int)floor(vi[1]+0.1); iy %= ny; iy = (iy < 0)? (iy+ny):iy;
      iu = (int)(abs(xeq[i][2]));

      gid[i] = POW2_IDX(ix, iy, iu);

      if (ix > xhi_loc)  xhi_loc = ix;
      if (ix < xlo_loc)  xlo_loc = ix;
      if (iy > yhi_loc)  yhi_loc = iy;
      if (iy < ylo_loc)  ylo_loc = iy;

      j++;
    }
  }

  nx_loc = xhi_loc-xlo_loc+1;
  ny_loc = yhi_loc-ylo_loc+1;

  // check whether we found all indices
  int *exists;
  natoms_loc = nx_loc*ny_loc*nu;
  exists = (int *) memory->smalloc(natoms_loc*sizeof(int),
				    "FixGFMDStatic::comp_map::exists");
  memset(exists, -1, natoms_loc*sizeof(int));
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int ix = IX_FROM_POW2_IDX(gid[i]);
      int iy = IY_FROM_POW2_IDX(gid[i]);
      int iu = IU_FROM_POW2_IDX(gid[i]);
      
      int j  = LOC_IDX(ix, iy, iu);

      if (exists[j] >= 0) {
	char errstr[1024];
	int k = exists[j];
	sprintf(errstr, "fix gfmd/static: comp_map: Index position "
		"%i x %i x %i alread defined for atom %i at position "
		"%f %f %f, while this atom %i is at position %f %f %f.",
		ix, iy, iu,
		k, xeq[k][0], xeq[k][1], xeq[k][2],
		i, xeq[i][0], xeq[i][1], xeq[i][2]);
	error->one(errstr);
      }
      exists[j] = i;
    }
  }
  for (int ix = 0; ix < nx_loc; ix++) {
    for (int iy = 0; iy < ny_loc; iy++) {
      for (int iu = 0; iu < nu; iu++) {
	int j  = LOC_IDX(ix, iy, iu);

	if (exists[j] < 0) {
	  char errstr[1024];
	  sprintf(errstr, "fix gfmd/static: comp_map: Atom with index "
		  "%i x %i x %i missing.", ix, iy, iu);
	  error->one(errstr);
	}
      }
    }
  }
  memory->sfree(exists);
  

  if (gfmdlog) {
    fprintf(gfmdlog,"\nMapping info computed from initial configuration.\n");
  }

  *out_nx = nx;
  *out_ny = ny;
}

/* ----------------------------------------------------------------------
 * private method, to read the Phi matrix from binary file.
 * Interpolation is done with bilinear for elements near the boundary
 * and bicubic elsewhere. 
 * --------------------------------------------------------------------*/

void FixGFMDStatic::readphi(char *phi_fn)
{
  error->all("readphi is currently broken. fix it");

#if 0
  int  Nx, Ny, Nucell, idim, ndim;
  double boltz, old2new = 1.;
  double svec_gfc[2][2], sb_gfc[ndof];
  FILE *gfc_in;
  int nread;

  fprintf(gfmdlog, "Reading stiffness kernel from '%s'.\n", phi_fn);

  gfc_in = fopen(phi_fn, "rb");
  if (gfc_in == NULL) error->all("fix gfmd/static: error while opening binary Phi file");

  //TASdebugout if (me == 0) fprintf(gfmdlog,"\nTAS For some reason, we're reading phi from file\n");

  nread = fread(&ndim,  sizeof(int),   1, gfc_in);
  if (nread != 1) error->one("fix gfmd/static: error while reading ndim from binary file");
  nread = fread(&Nx,    sizeof(int),   1, gfc_in);
  if (nread != 1) error->one("fix gfmd/static: error while reading Nx from binary file");
  nread = fread(&Ny,    sizeof(int),   1, gfc_in);
  if (nread != 1) error->one("fix gfmd/static: error while reading Ny from binary file");
  nread = fread(&Nucell,sizeof(int),   1, gfc_in);
  if (nread != 1) error->one("fix gfmd/static: error while reading Nucell from binary file");
  nread = fread(&boltz, sizeof(double),1, gfc_in);
  if (nread != 1) error->one("fix gfmd/static: error while reading boltz from binary file");
 
  if (ndim != sysdim){
    char str[MAXLINE];
    sprintf(str,"fix gfmd/static: System dimension from GFC is %d, current is %d", ndim, sysdim);
    if (me == 0){ fprintf(gfmdlog,"\n%s\n",str); fflush(gfmdlog); }
    error->all(str);
  }
  if (Nucell != nu){
    char str[MAXLINE];
    sprintf(str,"fix gfmd/static: # of atom per cell from GFC is %d, current is %d", Nucell, nu);
    if (me == 0){ fprintf(gfmdlog,"\n%s\n",str); fflush(gfmdlog); }
    error->all(str);
  }
  if (boltz != force->boltz){
    char str[MAXLINE];
    if (boltz == 1.){
      sprintf(str,"fix gfmd/static: Units used by GFC were LJ, conversion is not possible");
      if (me == 0){ fprintf(gfmdlog,"\n%s\n",str); fflush(gfmdlog); }
      error->all(str);
    } else {
      old2new = force->boltz / boltz;
      sprintf(str,"fix gfmd/static: Units used by GFC differ from current one, converted!");
      if (me == 0) fprintf(gfmdlog,"\n%s\n",str);
      error->warning(str);
    }
  }

  if (me == 0){
    fprintf(gfmdlog,"\nTo read  phi_q info from file: %s\n", phi_fn);
    fprintf(gfmdlog,"FFT mesh from GFC measurement: %d x %d\n", Nx, Ny);
    fprintf(gfmdlog,"FFT mesh from the current run: %d x %d\n", nx, ny);
  }

  if (nx == Nx && ny == Ny){ // Dimension from GFC run and current match, read in Phi directly
    double_complex cmdum;
    for (int idq=0; idq<xlo_loc*ny; idq++){ 
      for (idim=0; idim<ndof_sq; idim++){ // TAS read in all 9 components of Phi
        nread = fread(&cmdum, sizeof(double_complex), 1, gfc_in);
        if (nread != 1) error->one("fix gfmd/static: error while reading Phi from binary file");
      }
    }
    nread = fread(phi_q[0], sizeof(double_complex), mynq*ndof_sq, gfc_in); // TAS read in Phi(q) for all q's (although symmetric so only low qx's)
    if (nread != static_cast<int>(mynq*ndof_sq)) error->one("fix gfmd/static: error while reading Phi from binary file");
    if (me == 0) fprintf(gfmdlog,"\nphi_q read successfully from file: %s\n", phi_fn);

    if (xlo_loc==0 && mynq>0) phi_q0_rescale(phi_q[0]); // TAS phi_q[0] gets relabeled phi_q0 in phi_q0_rescale() 

  }else{ // Dimension from GFC run and current mismatch, interpolation needed!
    double_complex **Phi_in;  // read in phi_q from file
    Phi_in = create_2d_complex_array(memory, Nx*Ny, ndof_sq, "fix_gfmd:Phi_in");
    int idq = 0;
    for (int i=0; i<Nx; i++){
      for (int j=0; j<Ny; j++){
        for (idim=0; idim<ndof_sq; idim++){
          nread = fread(&Phi_in[idq][idim], sizeof(double_complex), 1, gfc_in);
          if (nread != 1) error->one("fix gfmd/static: error while reading Phi from binary file");
        }
        idq++;
      }
    }
    phi_q0_rescale(Phi_in[0]);  

    /* Now to interpolate the phi_q we needed!
       Bilinear interpolation is employed for boundary elements;
       while bicubic interpolation is used elsewhere. */
    int ix, iy, xP1, xM1, yP1, yM1, Idx, Ix, Iy;
    int ppx, ppy, pmx, pmy, mpx, mpy, mmx, mmy;
    double dx1, dx2, facx, facy;
    double_complex y[4], y1[4], y2[4], y12[4];
    double_complex *Phi_p1, *Phi_p2, *Phi_p12;
    Phi_p1  = new double_complex[(Nx+1)*(Ny+1)];
    Phi_p2  = new double_complex[(Nx+1)*(Ny+1)];
    Phi_p12 = new double_complex[(Nx+1)*(Ny+1)];
    dx1 = double(Nx)/2.;
    dx2 = double(Ny)/2.;
    
    for (idim=0; idim<ndof_sq; idim++){
      // get the gradients by finite element method
      for (Ix=0; Ix<=Nx; Ix++){
        for (Iy=0; Iy<=Ny; Iy++){
          
          int Cx = Ix%Nx, Cy = Iy%Ny;

          xP1 = (Cx+1)%Nx;
          xM1 = (Nx+Cx-1)%Nx;
          yP1 = (Cy+1)%Ny;
          yM1 = (Ny+Cy-1)%Ny;
          ppx = pmx = xP1;
          ppy = mpy = yP1;
          mpx = mmx = xM1;
          pmy = mmy = yM1;
          facx = 1.;
          facy = 1.;
         
          if (Ix == 0){
            xM1  = Cx;
            facx = 2.;
            mpx  = mmx = Cx;
            mpy  = mmy = Cy;
          } else if (Ix == Nx){
            xP1  = Cx;
            facx = 2.;
            ppx = pmx = Cx;
            ppy = pmy = Cy;
          }

          if (Iy == 0){
            yM1  = Cy;
            facy = 2.;
            pmx  = mmx = Cx;
            pmy  = mmy = Cy;
          } else if (Iy == Ny){
            yP1  = Cy;
            facy = 2.;
            ppx  = mpx = Cx;
            ppy  = mpy = Cy;
          }

          Idx = Ix * (Ny+1) + Iy;
          Phi_p1 [Idx] = (Phi_in[xP1*Ny+Cy ][idim] - Phi_in[xM1*Ny+Cy][idim]) * dx1 * facx;
          Phi_p2 [Idx] = (Phi_in[Cx*Ny+yP1 ][idim] - Phi_in[Cx*Ny+yM1][idim]) * dx2 * facy;
          Phi_p12[Idx] = (Phi_in[ppx*Ny+ppy][idim] - Phi_in[pmx*Ny+pmy][idim]
                       -  Phi_in[mpx*Ny+mpy][idim] + Phi_in[mmx*Ny+mmy][idim]) * dx1 * dx2 * facx * facy;
        }
      }

      // to do interpolation
      int idq = 0;
      for (ix=xlo_loc; ix<=xhi_loc; ix++){
        for (iy=0; iy<ny; iy++){
          Ix = (int)(double(ix)/double(nx)*Nx);
          Iy = (int)(double(iy)/double(ny)*Ny);
          xP1 = (Ix+1)%Nx;
          yP1 = (Iy+1)%Ny;

          y[0] = Phi_in[Ix*Ny+Iy][idim];
          y[1] = Phi_in[xP1*Ny+Iy][idim];
          y[2] = Phi_in[xP1*Ny+yP1][idim];
          y[3] = Phi_in[Ix*Ny+yP1][idim];

          xP1   = Ix+1;
          yP1   = Iy+1;
          y1[0] = Phi_p1[Ix *(Ny+1)+Iy ];
          y1[1] = Phi_p1[xP1*(Ny+1)+Iy ];
          y1[2] = Phi_p1[xP1*(Ny+1)+yP1];
          y1[3] = Phi_p1[Ix *(Ny+1)+yP1];

          y2[0] = Phi_p2[Ix *(Ny+1)+Iy];
          y2[1] = Phi_p2[xP1*(Ny+1)+Iy];
          y2[2] = Phi_p2[xP1*(Ny+1)+yP1];
          y2[3] = Phi_p2[Ix *(Ny+1)+yP1];

          y12[0] = Phi_p12[Ix *(Ny+1)+Iy];
          y12[1] = Phi_p12[xP1*(Ny+1)+Iy];
          y12[2] = Phi_p12[xP1*(Ny+1)+yP1];
          y12[3] = Phi_p12[Ix *(Ny+1)+yP1];

          bicuint(y, y1, y2, y12, double(Ix)/Nx, double(Ix+1)/Nx, double(Iy)/Ny,
                  double(Iy+1)/Ny, double(ix)/nx, double(iy)/ny, &phi_q[idq][idim]);
          idq++;
        }
      } // end of interpolation on current idim
    } // end of for (idim ...
     
    delete []Phi_p1;
    delete []Phi_p2;
    delete []Phi_p12;
    destroy_2d_complex_array(memory, Phi_in);
    if (me == 0) fprintf(gfmdlog,"\nphi_q interpolated from file: %s\n", phi_fn);
  }

  // to read surface lattice info from gfc measurements
  // in previous versions, this info is not stored.
  int info_read = 1;
  nread = fread(svec_gfc[0],sizeof(double),4,gfc_in);
  if (nread != 4){
    error->warning("fix gfmd/static: failed to read surface vector info from binary file");
    info_read = 0;
  }
  nread = fread(sb_gfc,sizeof(double),ndof, gfc_in);
  if (nread != static_cast<int>(ndof) ){
    error->warning("fix gfmd/static: failed to read surface basis info from binary file");
    info_read = 0;
  }
  fclose(gfc_in);

  // unit conversion for the elastic stiffness coefficients
  if (old2new != 1.){
    for (int idq=0; idq<mynq; idq++) {
      for (idim=0; idim<ndof_sq; idim++) {
	phi_q[idq][idim] *= old2new;
      }
    }
  }

  // compare equilibrium surface lattice based on gfc measurement and this run
  if (me == 0 && info_read){
    double sb_eq[ndof],d2o[sysdim];
    for (idim=0; idim<ndof; idim++) sb_eq[idim] = 0.;
    for (int ix=0; ix<nx; ix++){
      for (int iy=0; iy<ny; iy++){
        int idx = (ix*ny+iy)*nu;
        for (int iu=1; iu<nu; iu++){
          for (idim=0; idim<sysdim; idim++) d2o[idim] = xeq[idim][idx+iu] - xeq[idim][idx];
          domain->minimum_image(d2o);
          ndim = iu*sysdim;
          for (idim=0; idim<sysdim; idim++) sb_eq[ndim+idim] += d2o[idim];
        }
      }
    }
    for (idim=sysdim; idim<ndof; idim++) sb_eq[idim] /= nx*ny;

    fprintf(gfmdlog,"\nSurface vector from this run: [%lg %lg], [%lg %lg]\n",
      surfvec[0][0], surfvec[0][1], surfvec[1][0], surfvec[1][1]);
    fprintf(gfmdlog,"Surface vector from gfc  run: [%lg %lg], [%lg %lg]\n",
      svec_gfc[0][0], svec_gfc[0][1], svec_gfc[1][0], svec_gfc[1][1]);
    fprintf(gfmdlog,"Surface basis  from this run: ");
    for (idim=0; idim<ndof; idim++) fprintf(gfmdlog,"%lg ",sb_eq[idim]);
    fprintf(gfmdlog,"\nSurface basis  from gfc  run: ");
    for (idim=0; idim<ndof; idim++) fprintf(gfmdlog,"%lg ",sb_gfc[idim]);
    fprintf(gfmdlog, "\n");
  }
  // to reset xeq if required
  if (reset_xeq && info_read){
    double lx2now, ly2now, lx2gfc, ly2gfc, xs[3], dx2, dy2;
    lx2now = surfvec[0][0]*surfvec[0][0] + surfvec[0][1]*surfvec[0][1];
    ly2now = surfvec[1][0]*surfvec[1][0] + surfvec[1][1]*surfvec[1][1];
    lx2gfc = svec_gfc[0][0]*svec_gfc[0][0] + svec_gfc[0][1]*svec_gfc[0][1];
    ly2gfc = svec_gfc[1][0]*svec_gfc[1][0] + svec_gfc[1][1]*svec_gfc[1][1];
    dx2 = fabs(lx2now-lx2gfc); dy2 = fabs(ly2now-ly2gfc);
    if (dx2 > lx2now*2.5e-3 || dy2 > ly2now*2.5e-3){
      if (me == 0) error->warning("fix_gfmd: surface lattice from gfc and this run mismatch, xeq not reset!");
      return;
    }
    xs[0] = sqrt(lx2now/lx2gfc); xs[1] = sqrt(ly2now/ly2gfc);
    xs[2] = sqrt((lx2now+ly2now)/(lx2gfc+ly2gfc));
    for (idim=0; idim<ndof; idim++) sb_gfc[idim] *= xs[idim%sysdim];

    idx=0;
    for (int ix=0; ix<nx; ix++){
      for (int iy=0; iy<ny; iy++){
        ndim = 0;
        for (int iu=0; iu<nu; iu++){
          xeq[0][idx] = double(ix)*svec_gfc[0][0]+double(iy)*svec_gfc[1][0]+sb_gfc[ndim];
          xeq[1][idx] = double(iy)*svec_gfc[1][1]+sb_gfc[ndim+1];
          if (sysdim == 3) xeq[2][idx] = sb_gfc[ndim+2];
          idx++;
          ndim += sysdim;
        }
      }
    }
    if (me == 0) fprintf(gfmdlog,"\nEquilibrium positions reset based on surface lattice info from: %s\n",phi_fn);
  }
  return;
#endif
}



/* ----------------------------------------------------------------------
 * Private method, to set the acoustic frequencies to zero and reform
 * the Phi matrix at q = 0.
 * --------------------------------------------------------------------*/
// TAS: phi_q0_rescale() only called from function readphi()
void FixGFMDStatic::phi_q0_rescale(double_complex * phi_q0)
{
  double **phiq0, **vecq0, **inv_vecq0, **eigmat;
  double eigq0[ndof];
  int    ndim, idim, jdim;

  if (gfmdlog) {
    fprintf(gfmdlog,"\nTAS For some reason, we're rescaling phi_q0 from file\n");
  }

  if (ndof>sysdim){ // for system with basis, minimize acoustic eigenvalues

    phiq0     = memory->create_2d_double_array(ndof, ndof, "fix_gfmd:phiq0");
    vecq0     = memory->create_2d_double_array(ndof, ndof, "fix_gfmd:vecq0");
    inv_vecq0 = memory->create_2d_double_array(ndof, ndof, "fix_gfmd:inv_vecq0");
    eigmat    = memory->create_2d_double_array(ndof, ndof, "fix_gfmd:eigmat");
  
    ndim = 0;
    for (idim=0; idim<ndof; idim++){
      for (jdim=0; jdim<ndof; jdim++){
        phiq0[idim][jdim]  = creal(phi_q0[ndim++]);
        eigmat[idim][jdim] = 0.;
      }
    }
  
    if (jacobi(phiq0, ndof, eigq0, vecq0))
      error->one("fix gfmd/static: Max-iteration reached in getting eigenvalues of Phi(q=0)");
    if (gfmdlog){
      fprintf(gfmdlog, "\nEigenvalues of Phi(q0) from original gfc:");
      for (idim=0; idim<ndof; idim++) fprintf(gfmdlog," %lg", eigq0[idim]);
      fflush(gfmdlog);
    }
    for (idim=0; idim<sysdim; idim++) eigmat[idim][idim] = eigq0[idim]*1.e-10;
    for (idim=sysdim; idim<ndof; idim++) eigmat[idim][idim] = eigq0[idim];
    for (idim=0; idim<ndof; idim++){
      for (jdim=0; jdim<ndof; jdim++) inv_vecq0[idim][jdim] = vecq0[idim][jdim];
    }
  
    GaussJordan(ndof,inv_vecq0[0], error);
    MatMulMat(ndof,vecq0[0],eigmat[0],phiq0[0]);
    MatMulMat(ndof,phiq0[0],inv_vecq0[0],vecq0[0]); // Now vecq0 carries phiq0
  
    ndim = 0;
    for (idim=0; idim<ndof; idim++) {
      for (jdim=0; jdim<ndof; jdim++) {
	phi_q0[ndim++] = vecq0[idim][jdim];
      }
    }

    if (gfmdlog){
      jacobi(vecq0, ndof, eigq0, phiq0);
      fprintf(gfmdlog, "\nEigenvalues after rescaling original Phi:");
      for (idim=0; idim<ndof; idim++) fprintf(gfmdlog," %lg", eigq0[idim]);
      fprintf(gfmdlog,"\n");
    }

    memory->destroy_2d_double_array(phiq0);
    memory->destroy_2d_double_array(vecq0);
    memory->destroy_2d_double_array(inv_vecq0);
    memory->destroy_2d_double_array(eigmat);

  }else{ // for primitive systems, simply zero Phi(q=0).
    if (me == 0) fprintf(screen,"Entered phi_q0_rescale!!!\n");
    for (idim=0; idim<ndof_sq; idim++) {
      // TAS zero's all 9 Phi(q=0) elements
      phi_q0[idim] = 0.0;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
 * To output the elastic force, if keyword output is set.
 * --------------------------------------------------------------------*/

void FixGFMDStatic::end_of_step()
{
    /*
      FIXME!!!
  if (me == 0){
    char file_for[MAXLINE];
    FILE *fp;

    sprintf(file_for,"%s.%d", prefix, update->ntimestep);
    fp = fopen(file_for, "w");

    fprintf(fp,"# Elastic forces acting on GF layer, timestep=%d\n",update->ntimestep);
    fprintf(fp,"# Size Info: %lg %lg %lg %d %d\n", domain->xprd, domain->yprd,
      domain->xy, natoms, sysdim);
    if (sysdim == 3) fprintf(fp,"#index atom x y z fx fy fz\n");
    else fprintf(fp,"# index atom x y fx fy\n");

    for (int i=0; i<natoms; i++){
      idx  = static_cast<int>(UIrAll[i][sysdim]);
      itag = surf2tag[idx];
      fprintf(fp,"%d %d", idx, itag);
      for (int idim=0;idim<sysdim;idim++) fprintf(fp," %lg", UIrAll[i][idim]);
      for (int idim=0;idim<sysdim;idim++) fprintf(fp," %lg", FrAll[i][idim]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
    */
}

/* ----------------------------------------------------------------------
 * private method to compute the center-of-mass coordinates of GFMD atoms
 * at equilibrium positions
 * --------------------------------------------------------------------*/
 
void FixGFMDStatic::get_cm_from_xeq(double *xcm)
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

double FixGFMDStatic::memory_usage()
{
  double bytes = solver->memory_usage();

  // phi_q
  bytes += nxy_loc*ndof_sq*sizeof(double_complex);
  // u_xy
  bytes += nxy_loc*ndof*sizeof(double);
  // f_xy
  bytes += nxy*ndof*sizeof(double);

  return bytes;
}



