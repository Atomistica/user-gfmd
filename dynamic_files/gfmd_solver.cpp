#ifdef GFMD_CUFFT
#include <pthread.h>
#endif

#include <string.h>
#include <stdlib.h>

#include "pointers.h"

#include "comm.h"
#include "domain.h"
#include "linearalgebra.h"
#include "memory.h"
#include "mpi.h"

#include "surface_stiffness.h"

#include "gfmd_misc.h"
#include "gfmd_solver.h"

#include "gfmd_solver_static.h"
#include "gfmd_solver_dynamic.h"
#ifdef LMP_USER_CUDA 
#include "gfmd_solver_cuda.h"
#endif

using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))


/* ----------------------------------------------------------------------
 * Base class for 2d GFMD grids
 * --------------------------------------------------------------------*/

GFMDSolver::GFMDSolver(LAMMPS *lmp) : Pointers(lmp)
{

  strcpy(name, "none");

  u0 = NULL;
  nxwall = NULL;
  nywall = NULL;
  fxwall = NULL;
  fywall = NULL;
  /*
   * Get basic MPI information
   */
  me = comm->me;
  nprocs = comm->nprocs;
}


GFMDSolver::~GFMDSolver()
{
  if (u0)
    delete [] u0;
  if (nxwall){
    delete nxwall;
    delete fxwall;
  }
  if (nywall){
    delete nywall;
    delete fywall;
  }
}


void GFMDSolver::set_grid_size(int in_nx, int in_ny, int in_ndof)
{
  nx = in_nx;
  ny = in_ny;
  nu_ = in_ndof/3;
  ndof = in_ndof;

  ndof_sq = ndof*ndof;

  if (u0)
    delete u0;
  u0 = new double[ndof];

  /*
   * Check which part of the grid should be stored on this processor
   */

  xlo_loc = (int) round(nx*(domain->sublo[0]-domain->boxlo[0])/domain->xprd);
  xhi_loc = (int) round(nx*(domain->subhi[0]-domain->boxlo[0])/domain->xprd)-1;
  
  ylo_loc = (int) round(ny*(domain->sublo[1]-domain->boxlo[1])/domain->yprd);
  yhi_loc = (int) round(ny*(domain->subhi[1]-domain->boxlo[1])/domain->yprd)-1;
    
  if (!strcmp(this->get_name(), "dynamic")){
    nprocs_x = comm->procgrid[0];
    nprocs_y = comm->procgrid[1];
    nprocs_z = comm->procgrid[2];
    nwalls_x = nprocs_x+1;
    nwalls_y = nprocs_y+1;
    int lb_flag = get_lb_flag();
    int nmax = get_nmax();
    int nmin = get_nmin();

    int* procxwall = new int[nwalls_x];
    int* procywall = new int[nwalls_y];
    
    if (me==0) printf(" I see lb_flag = %i and nx = %i and ny = %i and nmax is %i and nmin is %i on %i procs with %i x walls and %i y walls\n",lb_flag, nx, ny, nmax, nmin, nprocs, nwalls_x, nwalls_y);
    if (lb_flag == 1 && nprocs > 4) {
      if (nprocs_z > 1) error->all(FLERR,"solver default/dynamic: Domain decomposition only works with 2D processor grid\n");
      
      if (me==0) {
	printf("Now attempting to do some load balancing\n");
	if (nxwall){
	  delete nxwall;
	  delete fxwall;
	}
	if (nywall){
	  delete nywall;
	  delete fywall;
	}
	nxwall = new int[nwalls_x];
	nywall = new int[nwalls_y];
	fxwall = new int[nwalls_x];
	fywall = new int[nwalls_y];
	
	//memory->create(nxwall, nwalls_x, "GFMDSolver::nxwall");
	//memory->create(nywall, nwalls_y, "GFMDSolver::nywall");
	
	memset(nxwall, 0, nwalls_x*sizeof(int));
	memset(nywall, 0, nwalls_y*sizeof(int));
	memset(fxwall, 0, nwalls_x*sizeof(int));
	memset(fywall, 0, nwalls_y*sizeof(int));
	
	nxwall[0] = 0;
	nywall[0] = 0;
	nxwall[nprocs_x] = nx-1;
	nywall[nprocs_y] = ny-1;
	
	int stride_x = nx/nprocs_x;
	int stride_y = ny/nprocs_y;
	for (int i=1; i<nprocs_x; i++){
	  nxwall[i] = i*stride_x;
	}
	for (int i=1; i<nprocs_y; i++){
	  nywall[i] = i*stride_y;
	}

	ntot = layer_sum(0,nx-1,0,ny-1,nx,ny,nmax,nmin);
	printf("nlayers = %i | ", int(ntot));
	nxtarget = float(ntot)/nprocs_x;
	nytarget = float(ntot)/nprocs_y;

	printf("Per processor target = %f \n",float(ntot)/nprocs);
	printf("nxtarget = %f nytarget = %f\n\n",nxtarget, nytarget);
	
	for (int i=0; i<nwalls_x; i++){
	  printf("xi = %i wall at %i | ", i, nxwall[i]);
	}
	printf("\n");
	for (int i=0; i<nwalls_y; i++){
	  printf("yi = %i wall at %i | ", i, nywall[i]);
	}
	printf("\n\n");
	
	int nsum, nx1, ny1, fxsum, fysum;
	double fx1, fy1;
	for (int count=0; count < 4; count++){
	  // Collect layer info
	  for (int i=0; i<nprocs_x; i++){
	    //printf("\n count = %i | i = %i | nxwall = %i %i | ", count, i, nxwall[i],nxwall[i+1]);
	    nx1 = layer_sum(nxwall[i], nxwall[i+1]-1, 0, ny-1, nx, ny, nmax, nmin);
	    fx1 = int((nxtarget-nx1)/nxtarget*10);
	    if (SIGN(fx1)*fx1>1) fx1 = SIGN(fx1);
	    if (i>0) fxwall[i] -= int(fx1);
	    if (i<nprocs_x) fxwall[i+1] += int(fx1);
	    //printf("layers = %i | fx = %f \n", nx1, fx1);
	  }
	  for (int i=0; i<nprocs_y; i++){
	    //printf("\n count = %i | i = %i | nywall = %i %i | ", count, i, nywall[i],nywall[i+1]);
	    ny1 = layer_sum(0,nx-1,nywall[i], nywall[i+1]-1, nx, ny, nmax, nmin);
	    fy1 = int((nytarget-ny1)/nytarget*10);
	    if (SIGN(fy1)*fy1>1) fy1 = SIGN(fy1);
	    if (i>0) fywall[i] -= int(fy1);
	    if (i<nprocs_y) fywall[i+1] += int(fy1);
	    //printf("layers = %i | fy = %f \n", ny1, fy1);
	  }
	  
	  // Do wall integration
	  fxsum = 0;
	  int step;
	  for (int i=1; i<nprocs_x; i++){
	    step = SIGN(fxwall[i]);
	    nxwall[i] += step;
	    fxsum += step;
	    fxwall[i] = 0;
	    if (nxwall[i]<=nxwall[i-1]) nxwall[i] = nxwall[i-1]+1;
	    if (nxwall[i]>=nxwall[i+1]) nxwall[i] = nxwall[i+1]-1;
	  }
	  fysum = 0;
	  for (int i=1; i<nprocs_y; i++){
	    step = SIGN(fywall[i]);
	    nywall[i] += step;
	    fysum += step;
	    fywall[i] = 0;
	    if (nywall[i]<=nywall[i-1]) nywall[i] = nywall[i-1]+1;
	    if (nywall[i]>=nywall[i+1]) nywall[i] = nywall[i+1]-1;
	  }
	  
	  // Check output
	  /*
	  for (int i=0; i<nwalls_x; i++){
	    printf("xi = %i wall at %i | ", i, nxwall[i]);
	  }
	  printf("and fxsum = %i \n", fxsum);
	  for (int i=0; i<nwalls_y; i++){
	    printf("yi = %i wall at %i | ", i, nywall[i]);
	  }
	  printf("and fysum = %i \n", fysum);
	  */
	    
	  int layer_count = 0;
	  printf("\ncount = %i\n",count);
	  for (int i=0;i<nprocs_x;i++){
	    for (int j=0;j<nprocs_y;j++){
	      layer_count = layer_sum(nxwall[i],nxwall[i+1]-1,nywall[j], nywall[j+1]-1, nx, ny, nmax, nmin);
	      printf("x: %i %i | y: %i %i | n = %i\n", nxwall[i], nxwall[i+1]-1, nywall[j], nywall[j+1]-1, layer_count);
	      //printf("%i\n", layer_count);
	    }
	  }
	  printf("\n");
	  
	  // Pack Buffer
	  for (int i=0;i<nwalls_x;i++){
	    procxwall[i] = nxwall[i];
	  }
	  for (int i=0;i<nwalls_y;i++){
            procywall[i] = nywall[i];
          }

	}

      }
      MPI_Bcast(procxwall, nwalls_x, MPI_INT, 0, world);
      MPI_Bcast(procywall, nwalls_y, MPI_INT, 0, world);

      int rank_x = me % nprocs_x;
      int rank_y = me / nprocs_x;
      xlo_loc_ks = procxwall[rank_x];
      xhi_loc_ks = procxwall[rank_x+1]-1;
      ylo_loc_ks = procywall[rank_y];
      yhi_loc_ks = procywall[rank_y+1]-1;
      printf("I am %i with bounds x: %i %i | y: %i %i\n\n", me, xlo_loc_ks, xhi_loc_ks, ylo_loc_ks, yhi_loc_ks);
    }
    else {
      xlo_loc_ks = xlo_loc;
      xhi_loc_ks = xhi_loc;

      ylo_loc_ks = ylo_loc;
      yhi_loc_ks = yhi_loc;
      
    }

    nx_loc_ks  = xhi_loc_ks-xlo_loc_ks+1;
    ny_loc_ks  = yhi_loc_ks-ylo_loc_ks+1;
    
    nxy_loc_ks = nx_loc_ks*ny_loc_ks;
    
  }

  int nxlo0 = 0, nylo0 = 0;

  if (xlo_loc == 0) nxlo0 = 1;
  if (ylo_loc == 0) nylo0 = 1;

  int nxlo0_red, nylo0_red;
  MPI_Reduce(&nxlo0, &nxlo0_red, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nylo0, &nylo0_red, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (me == 0) {
    if (nxlo0_red != comm->procgrid[1]) {
      char errstr[1024];
      sprintf(errstr, "[GFMDSolver::GFMDSolver] %i processors have x-grid "
	      "index 0, while we have a decomposition of %i processors in "
	      "y-direction.", nxlo0_red, comm->procgrid[1]);
      error->one(FLERR,errstr);
    }
    if (nylo0_red != comm->procgrid[0]) {
      char errstr[1024];
      sprintf(errstr, "[GFMDSolver::GFMDSolver] %i processors have y-grid "
	      "index 0, while we have a decomposition of %i processors in "
	      "x-direction.", nylo0_red, comm->procgrid[0]);
      error->one(FLERR,errstr);
    }
  }

  nx_loc  = xhi_loc-xlo_loc+1;
  ny_loc  = yhi_loc-ylo_loc+1;

  nxy_loc = nx_loc*ny_loc;

  gammai_ = -1;
  if (xlo_loc <= 0 && 0 <= xhi_loc && ylo_loc <= 0 && 0 <= yhi_loc) {
    gammai_ = -xlo_loc*ny_loc - ylo_loc;
  }
}

int GFMDSolver::layer_sum(int xlo_loc, int xhi_loc, int ylo_loc,
			 int yhi_loc, int nx, int ny, int nmax, int nmin)
{
  int nsum = 0;
  int ntmp;
  int max_dim =  (nx>=ny) ? (nx) : (ny);

  //printf(" xlo_loc = %i | xhi_loc = %i | ylo_loc = %i | yhi_loc = %i | nx = %i | ny = %i \n",xlo_loc,xhi_loc, ylo_loc, yhi_loc, nx,ny);
  //if (xhi_loc == nx-1) xhi_loc += 1;
  //if (yhi_loc == ny-1) yhi_loc += 1;
  double qmag, qx, qy;
  for (int i = xlo_loc; i <= xhi_loc; i++) {
    qx = (i <= int(nx/2)) ?
      (2.0*M_PI*i/nx) : (2.0*M_PI*(i-nx)/nx);
    for (int j = ylo_loc; j <= yhi_loc; j++) {
      qy = (j <= int(ny/2)) ?
        (2.0*M_PI*j/ny) : (2.0*M_PI*(j-ny)/ny);
      
      qmag = sqrt(qx*qx+qy*qy);
      if( qmag < 0.000001) {
	ntmp = nmax;
      }
      else {
        ntmp = round(nmax/(qmag/(2.0*M_PI)*max_dim));
        if (ntmp<nmin) ntmp = nmin;
      }
      //printf("xx %i %f %i\n",comm->me, qmag, ntmp);
      nsum += ntmp;
    }
  }

  /*
  int nsum = 0;
  int ntmp = 0;

  for (int i = xlo_loc; i < xhi_loc; i++) {
    double qx = (i <= int((nx)/2)) ?
      (1.0*(i)/nx) : (1.0*(i-nx)/nx);
    for (int j = ylo_loc; j < yhi_loc; j++) {
      double qy = (j <= int((ny)/2)) ?
	(1.0*(j)/ny) : (1.0*(j-ny)/ny);
      
      double qmag = sqrt(qx*qx+qy*qy);
      if ( qmag < 0.000001){
	ntmp = nmax;
      }
      else {
	ntmp = round(nmax/qmag);
          }
      if (ntmp<nmin) ntmp = nmin;
      nsum += ntmp;
    }
    }*/
    return nsum;
}

double **GFMDSolver::create_double_buffer(const char *name)
{
  double **tmp2d;
  return memory->create(tmp2d, MAX(1, nxy_loc), ndof, name);
}


void GFMDSolver::destroy_double_buffer(double **buffer)
{
  memory->destroy(buffer);
}


double_complex **GFMDSolver::create_complex_buffer(const char *name)
{
  double_complex **tmp2d;
  return memory->create(tmp2d, MAX(1, nxy_loc), ndof, name);
}


void GFMDSolver::destroy_complex_buffer(double_complex **buffer)
{
  memory->destroy(buffer);
}


double_complex **GFMDSolver::create_complex_operator_buffer(const char *name)
{
  double_complex **tmp2d;
  return memory->create(tmp2d, MAX(1, nxy_loc), ndof_sq, name);
}


void GFMDSolver::destroy_complex_operator_buffer(double_complex **buffer)
{
  memory->destroy(buffer);
}


double GFMDSolver::memory_usage()
{
  return 0.0;
}


void GFMDSolver::dump_stiffness()
{
  char errstr[1024];
  sprintf(errstr, "GFMDSolver::dump_stiffness: Solver '%s' does not support "
	  "dumping of stiffness coefficients.", name);
  error->all(FLERR,errstr);
}


void GFMDSolver::dump_greens_function()
{
  char errstr[1024];
  sprintf(errstr, "GFMDSolver::dump_greens_function: Solver '%s' does not "
	  "support dumping of Green's functions coefficients.", name);
  error->all(FLERR,errstr);
}



namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Instantiate a GFMD solver according to keyword arguments
 * --------------------------------------------------------------------*/

GFMDSolver *gfmd_solver_factory(char *keyword, LAMMPS *lmp, int narg,
				int *carg, char **arg)
{
  GFMDSolver *solver = NULL;

  char solver_name[128];
  if (lmp->suffix_enable) {
    if (keyword) {
      snprintf(solver_name, 128, "%s/%s", keyword, lmp->suffix);
    }
    else {
      snprintf(solver_name, 128, "static/%s", lmp->suffix);
    }
  }
  else {
    if (keyword) {
      strncpy(solver_name, keyword, 128);
    } 
    else {
      strcpy(solver_name, "static");
    }
  }

  if (!strcmp(solver_name, "static")) {
    solver = new GFMDSolverStatic(lmp, narg, carg, arg);
  }
  else if (!strcmp(solver_name, "dynamic")) {
    solver = new GFMDSolverDynamic(lmp, narg, carg, arg);
  }
  
#ifdef LMP_USER_CUDA
  else if (!strcmp(solver_name, "static/cuda")) {
    solver = new GFMDSolverStaticCUDA(lmp);
  }
#endif
  else {
    lmp->error->one(FLERR,"Unknown solver name encountered.");
  }

  if (solver) {
    if (strcmp(solver->get_name(), solver_name)) {
      char errstr[120];

      sprintf(errstr, "gfmd_solver_factory: Internal error: solver_name '%s' "
	      "and solver class '%s' name mismatch.", solver_name,
	      solver->get_name());
      lmp->error->all(FLERR,errstr);
    }
  }

  return solver;
}

}
