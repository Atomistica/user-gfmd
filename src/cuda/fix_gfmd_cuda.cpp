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
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>

#include "mpi.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "fix_gfmd_cuda.h"
#include "force.h"
#include "group.h"
#include "lattice.h"
#include "linearalgebra.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"
#include "cuda.h" // note cuda_data.h etc included in cuda.h
#include "cuda_modify_flags.h"

#include "gfmd_grid.h"
#include "gfmd_solver.h"
#include "surface_stiffness.h"

#include "fix_gfmd_cuda_cu.h"

#ifdef GFMD_PROFILING
#include "cuda_profiler_api.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

#include "USER-GFMD/REV"

#define MAXLINE   256

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

#define nearbyint(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define rint(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define round(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))
#define lround(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))

#define modulo(x, y)  ( x < 0 ? x-((x-y+1)/y)*y : x-(x/y)*y )
#define fmodulo(x, y)  ( x < 0 ? x-floor((x-y+1)/y)*y : x-floor(x/y)*y )

/* ---------------------------------------------------------------------- */

FixGFMDCuda::FixGFMDCuda(LAMMPS *lmp, int narg, char **arg)
  : FixGFMD(lmp, narg, arg)
{
  cuda = lmp->cuda;

  if (nprocs > 1)
    error->all(FLERR,"GFMD and CUDA work on single CPU/GPU combinations only "
	       "for now.");
  cu_u_xy = NULL;
  cu_f_xy = NULL;
  cu_force_is_known = NULL;
} // end of constructor

/* ---------------------------------------------------------------------- */

FixGFMDCuda::~FixGFMDCuda()
{
  // If pointer is NULL, flags that memory is not available
  delete cu_u_xy; cu_u_xy=NULL;
  delete cu_f_xy; cu_f_xy=NULL;

  delete cu_force_is_known; cu_force_is_known=NULL;
}

/* ---------------------------------------------------------------------- */

int FixGFMDCuda::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE_CUDA;
  mask |= POST_FORCE_CUDA;
  mask |= THERMO_ENERGY_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGFMDCuda::init()
{
  FixGFMD::init();

  // displacements
  if (cu_u_xy == NULL) cu_u_xy = new cCudaData<double, X_FLOAT, xy>((double*)u_xy, rnu*ndof/nu,
					       nxy_loc, 0, true);
  // forces
  if (cu_f_xy == NULL) cu_f_xy = new cCudaData<double, X_FLOAT, xy>((double*)f_xy, rnu*ndof/nu,
					       nxy_loc, 0, true);

  // auxiliary buffer
  if (cu_force_is_known == NULL) cu_force_is_known = new cCudaData<bool, bool, x>(force_is_known,
						   atom->nmax, 0, 0, true);
}

/* ---------------------------------------------------------------------- */

void FixGFMDCuda::setup(int vflag)
{
  FixGFMD::setup(vflag);

  // need to download possibly changed arrays because VerletCuda::run will
  // upload again
  cuda->cu_f->download();
  cuda->cu_gid->download();
}

/* ---------------------------------------------------------------------- */

void FixGFMDCuda::pre_force(int vflag)
{
#ifdef GFMD_PROFILING
  cudaProfilerStart();
#endif

  X_FLOAT *dev_x = static_cast<X_FLOAT*>(cuda->cu_x->dev_data());
  X_FLOAT *dev_xeq = static_cast<X_FLOAT*>(cuda->cu_xeq->dev_data());
  //printf("ORIG devloc of x\n"); cuda->cu_x->printaddr(); // TAS
  
#ifdef DEBUG_TAS
  printf("cu_u_xy:"); 
  CudaData_printmyaddress(cu_u_xy->dev_data());// TAS  cu_u_xy->printaddr();
  printf("cu_f_xy:"); 
  CudaData_printmyaddress(cu_f_xy->dev_data());// TAS
  fflush(stdout);
#endif DEBUG_TAS

  int *dev_mask = static_cast<int*>(cuda->cu_mask->dev_data());
  int *dev_gid = static_cast<int*>(cuda->cu_gid->dev_data());

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double xprd_half = domain->xprd_half;
  double yprd_half = domain->yprd_half;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int natoms_cur;

  /* --- */

  double previous_time = MPI_Wtime();

  if (nmax < atom->nmax) {
    nmax = atom->nmax;
    force_is_known = (bool *) 
      memory->srealloc(force_is_known, nmax*sizeof(bool),
                       "FixGFMDCuda:force_is_known");
    delete cu_force_is_known;
    // auxiliary buffer
    cu_force_is_known = new cCudaData<bool, bool, x>(force_is_known, atom->nmax,
						     0, 0, true);
  }

  double dxcm[3];
  xcm(igroup, masstotal, dxcm);

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

  epot = 0.0;
  epot_flag = false;

  /*
   * get U(r) on local proc
   */
  ngfmd_loc = 0;
  natoms_cur = 0;
  if (domain->triclinic == 0) {
    // for orthogonal lattice
    if (fix_gfmd_cuda_list_to_grid(xprd, yprd, xprd_half, yprd_half, groupbit,
				   nlocal, nall, nmax, dev_mask, dev_x, dev_gid,
				   dev_xeq, dxshift, dyshift, nx, ny, xlo_loc,
				   ylo_loc, nx_loc, ny_loc, nxy_loc,
				   static_cast<double*>(cu_u_xy->dev_data()),
				   natoms_cur, ngfmd_loc) != cudaSuccess) {
      char errstr[1024];

      sprintf(errstr, "fix_gfmd_cuda_list_to_grid failed: %s.",
	      cudaGetErrorString(cudaGetLastError()));
      error->one(FLERR,errstr);
    }
  }
  else { 
    // for non-orthogonal lattice
    error->all(FLERR,"The non-orthogonal GFMD case is not implemented. "
	       "Fix it!");
  }  // TAS end of orthogonal/ non statements 

  if (natoms_cur != nxy_loc*rnu) {
    char errstr[1024];
    sprintf(errstr, "fix gfmd: natoms_cur != natoms_loc "
            "(%i != %i)", natoms_cur, nxy_loc*rnu);
  }
#ifdef DEBUG_TAS
  printf("cu_u_xy:"); 
  CudaData_printmyaddress(cu_u_xy->dev_data());// TAS  cu_u_xy->printaddr();
  printf("cu_f_xy:"); 
  CudaData_printmyaddress(cu_f_xy->dev_data());// TAS
  fflush(stdout);
#endif

  xshift = cur_xshift;
  yshift = cur_yshift;

  solver->pre_force(cu_u_xy->dev_data(), cu_f_xy->dev_data());

  pre_time += MPI_Wtime() - previous_time;

  pre_flag = true;

#ifdef GFMD_PROFILING
  cudaProfilerStop();
#endif
}


/* ---------------------------------------------------------------------- */

void FixGFMDCuda::post_force(int vflag)
{
#ifdef GFMD_PROFILING
  cudaProfilerStart();
#endif

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  /* --- */

  double previous_time = MPI_Wtime();

  if (!pre_flag)
    error->one(FLERR,"FixGFMDCuda: pre_force needs to be called before "
               "post_force");

  if (dumpq_every > 0 && update->ntimestep % dumpq_every == 0) {
    char dump_prefix[1024];
    sprintf(dump_prefix, "%s.%li", prefix, update->ntimestep);
    epot += solver->post_force(cu_u_xy->dev_data(), cu_f_xy->dev_data(),
    			       dump_prefix);
  }
  else 
    epot += solver->post_force(cu_u_xy->dev_data(), cu_f_xy->dev_data(), NULL);

  /*
   * Copy f_xy from grid back to force list
   */
  grid_to_list_cuda();
  fsum_flag = false;

  /*
   * Dump displacements, forces and energies
   */
  if (dumpr_every > 0 && update->ntimestep % dumpr_every == 0) {
    cu_u_xy->download(false);
    cu_f_xy->download(false);

#if 0
    for (int j = 0; j < ndof; j++) {
      bool differs = false;
      for (int i = 0; i < nxy_loc; i++) {
	//	printf("%e %e %e\n", f_xy[j][i]/256, u_xy[j][i], abs(f_xy[j][i]*256-u_xy[j][i]));
	if (abs(f_xy[j][i]/256-u_xy[j][i]) > 1e-12)
	  //if (abs(f_xy[j][i]) > 1e-12)
	  differs = true;
      }
      if (differs) printf("dim = %i differs\n", j);
    }
#endif
    
    dump();
  }

  /*
   * Remove force on the center of mass
   */
  if (fix_center[0] || fix_center[1] || fix_center[2]) {
				error->all(FLERR,"fix_center not (yet) supported.");
				/*
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
				*/
  }

  /*
   * Gather total energy on all processors
   */
  double tmp;
  MPI_Allreduce(&epot, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  epot = tmp; //if (epot > 10000) {printf("caught a big KE\n"); exit(999);}// TAS

  epot_flag = true;
  pre_flag  = false;

  post_time += MPI_Wtime() - previous_time;

#ifdef GFMD_PROFILING
  cudaProfilerStop();
#endif
} // end post_force(int)

/* ----------------------------------------------------------------------
 * copy force data from the grid f_xy back to the f_i list that stores
 * the forces per atom. f_i also contains forces on ghost atoms while
 * f_xy exclusively stores information on grid atoms.
 * --------------------------------------------------------------------*/

void FixGFMDCuda::grid_to_list_cuda()
{
#if 0
  int *gid   = atom->gid;
  int *mask  = atom->mask;
#else
  int *mask  = static_cast<int*>(cuda->cu_mask->dev_data());
  int *gid   = static_cast<int*>(cuda->cu_gid->dev_data());
#endif
  int nlocal = atom->nlocal;
  int nall   = nlocal + atom->nghost;

  F_FLOAT *dev_f = static_cast<F_FLOAT*>(cuda->cu_f->dev_data());
  double *dev_f_xy = static_cast<double*>(cu_f_xy->dev_data());
  bool *dev_force_is_known = static_cast<bool*>(cu_force_is_known->dev_data());

  int natoms_cur;

  /*
   * add elastic force to local atoms
   */
  natoms_cur = 0;
  fsum_loc[0] = fsum_loc[1] = fsum_loc[2] = 0.0;
  for (int i = 0; i < nall; i++) {
    force_is_known[i] = false;
  }

  // dev_f is the force list; dev_f_xy is the force grid
  if (fix_gfmd_cuda_grid_to_list(groupbit, nlocal, nall, nmax, mask, dev_f,
				 gid, nx_loc, ny_loc, nxy_loc, xlo_loc, ylo_loc,
				 xhi_loc, yhi_loc, dev_f_xy, dev_force_is_known,
				 natoms_cur, fsum_loc)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "fix_gfmd_cuda_grid_to_list failed: %s.",
	    cudaGetErrorString(cudaGetLastError()));
    error->one(FLERR,errstr);
  }

#if 0
  //  cu_f_xy->download(true);
  //  cuda->cu_f->download(true);

  /*
   * communicate ghost data back to this processor
   */
  natoms_loc_recv = 0;
  comm->reverse_comm_fix(this);
  natoms_cur += natoms_loc_recv;

  cu_f_xy->upload(true);

  /*
   * Check if anything went wrong
   */
  if (natoms_cur != ngfmd_loc) {
    double **x = atom->x;
    double **xeq = atom->xeq;
    int *tag = atom->tag;

    gid = atom->gid;
    mask = atom->mask;

    cuda->cu_x->download(false);
    cuda->cu_xeq->download(false);
    cuda->cu_tag->download(false);
    cuda->cu_gid->download(false);
    cuda->cu_mask->download(false);

    double xcm[3];
    group->xcm(igroup, masstotal, xcm);

    double *u0 = solver->get_u0();

    char loststr[1024];
    strcpy(loststr, "");
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (!force_is_known[i]) {
          char newstr[1024];

          int j = gid[i];

          int ix = IX_FROM_POW2_IDX(j);
          int iy = IY_FROM_POW2_IDX(j);
          int iu = IU_FROM_POW2_IDX(j);

          sprintf(newstr, "%i (%i %i %i - %f %f %f), ", tag[i],
                  ix, iy, iu, x[i][0]-xeq[i][0], x[i][1]-xeq[i][1],
                  x[i][2]-xeq[i][2]);
          if (strlen(loststr) + strlen(newstr) < 1024)
            strcat(loststr, newstr);
        }
      }
    }

    char errstr[2048];
    if (strlen(loststr) >= 2) {
      loststr[strlen(loststr)-2] = 0;

      sprintf(errstr, "fix gfmd: natoms_cur != ngfmd_loc (%i != %i). Are "
              "the initial grid positions set correctly in the input file? "
              "Overall GFMD layer center of mass is %f %f %f. q=0 "
              "displacement is %f %f %f. The atoms that were lost are %s.",
              natoms_cur, ngfmd_loc, xcm[0], xcm[1], xcm[2], u0[0], u0[1],
              u0[2], loststr);
    }
    else {
      sprintf(errstr, "fix gfmd: natoms_cur != ngfmd_loc (%i != %i), but "
              "apparently no atoms are lost. Weird.", natoms_cur, ngfmd_loc);
    }
    error->one(FLERR,errstr);
  }
#endif
}

/* ----------------------------------------------------------------------
 * compute the center-of-mass coordinates
 * --------------------------------------------------------------------*/

void FixGFMDCuda::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = group->bitmask[igroup];

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double *x = static_cast<double*>(cuda->cu_x->dev_data());
  int *mask = static_cast<int*>(cuda->cu_mask->dev_data());
  int *type = static_cast<int*>(cuda->cu_type->dev_data());
  double *mass = static_cast<double*>(cuda->cu_mass->dev_data());
  double *rmass = NULL;
  if (cuda->cu_rmass) {
    rmass = static_cast<double*>(cuda->cu_rmass->dev_data());
  }
  tagint *image = static_cast<tagint*>(cuda->cu_image->dev_data());
  int nlocal = atom->nlocal;

  // Have each MPI process find its contribution to x,y,z center-of-mass
  double cmone[3];  
  if (fix_gfmd_cuda_xcm(groupbit, xprd, yprd, zprd, nlocal, nmax,
			mask, type, mass, rmass, x, image, cmone)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "fix_gfmd_cuda_xcm failed: %s.",
	    cudaGetErrorString(cudaGetLastError()));
    error->one(FLERR,errstr);
  }

  // Reduce calculation of center-of-mass
  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }

}
