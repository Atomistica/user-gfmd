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
