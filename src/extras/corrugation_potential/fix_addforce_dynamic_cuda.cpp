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
#include "fix_addforce_dynamic.h" // because i want to inherit from here
#include "fix_addforce_dynamic_cuda.h"
#include "fix_addforce_dynamic_cuda_cu.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "variable.h"
#include "error.h"
#include "cuda_modify_flags.h"
#include "cuda.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;


// TAS developing by comparing to fix_viscous_cuda.cpp
/* ---------------------------------------------------------------------- */

FixAddForceDynamicCuda::FixAddForceDynamicCuda(LAMMPS *lmp, int narg, char **arg) :
  FixAddForceDynamic(lmp, narg, arg)
{
  cuda = lmp->cuda;
  if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");
}

/* ---------------------------------------------------------------------- */

FixAddForceDynamicCuda::~FixAddForceDynamicCuda()
{
}

/* ---------------------------------------------------------------------- */

int FixAddForceDynamicCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */
//
//void FixAddForceDynamicCuda::min_setup(int vflag)
//{
//  Cuda_FixAddForceDynamicCuda_Init(&cuda->shared_data);
//  post_force(vflag);
//}
//
///* ---------------------------------------------------------------------- */
//
void FixAddForceDynamicCuda::setup(int vflag)
{
//  // Please note: atom->ntypes+1 was changed to 1.  Seems like maybe should make 1d instead of 2d
//  // Also, maybe the fact that gamma is an arbitary-sized vector means that they used a cu_ variable instead of something simpler
//   if(not cu_reflxcenter_sigma)
//   cu_reflxcenter_sigma = new cCudaData<double, F_FLOAT, x> (reflxcenter_sigma,1);
//   if(not cu_reflycenter_sigma)
//   cu_reflycenter_sigma = new cCudaData<double, F_FLOAT, x> (reflycenter_sigma,1);
//   if(not cu_radspersigma_x)
//   cu_radspersigma_x = new cCudaData<double, F_FLOAT, x> (radspersigma_x,1);
//   if(not cu_radspersigma_y)
//   cu_radspersigma_y = new cCudaData<double, F_FLOAT, x> (radspersigma_y,1);
//   if(not cu_aradsq)
//   cu_aradsq = new cCudaData<double, F_FLOAT, x> (aradsq,1);
//   if(not cu_border_dist)
//   cu_border_dist = new cCudaData<double, F_FLOAT, x> (border_dist,1);
//   if(not cu_r_overEstar)
//   cu_r_overEstar = new cCudaData<double, F_FLOAT, x> (r_overEstar,1);
//

  //if(not cu_foriginaladf)
  //cu_foriginaladf = new cCudaData<double, F_FLOAT, x> (foriginaladf,1);

  //Cuda_FixAddForceDynamicCuda_Init(&cuda->shared_data);

   //   cu_cu_foriginaladf->upload();

   //   cu_reflxcenter_sigma->upload();
   //   cu_reflycenter_sigma->upload();
   //   cu_radspersigma_x->upload();
   //   cu_radspersigma_y->upload();
   //   cu_aradsq->upload();
   //   cu_border_dist->upload();
   //   cu_cu_r_overEstar->upload();


  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    if (logfile) fprintf(logfile, " adf not defined for respa.  exit. \n");
    printf(" adf not defined for respa.  exit. \n");
    exit(999);
  }
  // need to download possibly changed arrays because VerletCuda::run will
  // upload again
  cuda->cu_f->download();

}

/* ---------------------------------------------------------------------- */

void FixAddForceDynamicCuda::post_force(int vflag)
{

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double xcenter_sigma = xprd / 2;
  double ycenter_sigma = yprd / 2;  
  double acontact_sigma = nndist_ * acontact_nn_;
  
  double truexcenter_sigma = xcenter_sigma + (shiftx_nn_ * nndist_);
  double reflxcenter_sigma = truexcenter_sigma - xprd*floor(truexcenter_sigma / xprd);

  double trueycenter_sigma = ycenter_sigma + (shifty_nn_ * nndist_);
  double reflycenter_sigma = trueycenter_sigma - yprd*floor(trueycenter_sigma / yprd);

  // Define contact potential
  double radspersigma_x = 2.0*M_PI / nndist_;
  double radspersigma_y = 2.0*M_PI / nndist_;
  double border_dist = 0.5 * nndist_;

  // To device's array of atoms, add a force that varies in space
  Cuda_FixAddForce_DynamicCuda_PostForce(&cuda->shared_data, groupbit,
					reflxcenter_sigma, 
					reflycenter_sigma,
					radspersigma_x,
					radspersigma_y,
		    			acontact_sigma, 
		    			border_dist,
		    			r_overEstar);
  // Rather than passing each post_force, could assign to shared memory at init()


}
