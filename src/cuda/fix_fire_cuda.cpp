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

/* ----------------------------------------------------------------------
   Contributing author: Lars Pastewka (JHU/Fh-IWM)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <cuda_runtime.h>

#include "atom.h"
#include "cuda.h"
#include "cuda_modify_flags.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "respa.h"

#include "fix_fire_cuda.h"
#include "fix_fire_cuda_cu.h"
#include "fix_nve_cuda_cu.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

#define CU_CHECK(x)							\
  {									\
  cudaError err = x;							\
  if (err != cudaSuccess) {						\
    char errstr[1024];							\
    sprintf(errstr, "CUDA error: %s.",					\
	    cudaGetErrorString(err));					\
    error->one(FLERR,errstr);						\
  }									\
  }


/* ---------------------------------------------------------------------- */

FixFIRECuda::FixFIRECuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;

  if(cuda == NULL)
    error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' "
	       "acceleration. Provide '-c on' as command-line argument to "
	       "LAMMPS..");

  minsteps = 10;
  incfac = 1.2;
  decfac = 0.5;
  mix_in = 0.1;
  mixdec = 0.99;
  max_dt = 0.05;
  xlimit = -1; // Not implemented in the cuda version of fix_fire
  delayuploadfrequency = 5; // Not implemented

  if (narg == 10) {
    minsteps = force->inumeric(FLERR,arg[3]);
    incfac = force->numeric(FLERR,arg[4]);
    decfac = force->numeric(FLERR,arg[5]);
    mix_in = force->numeric(FLERR,arg[6]);
    mixdec = force->numeric(FLERR,arg[7]);
    max_dt = force->numeric(FLERR,arg[8]);
    xlimit = force->numeric(FLERR,arg[9]); // Not implemented in the cuda version of fix_fire
  }
  else {
    if (narg != 3) error->all(FLERR,"Illegal fix fire command. Parameters are: "
			      "fix ID group-ID fire [minsteps incfac decfac "
			      "min_in mixdec max_dt xlimit(NOTIMPLEMENTED)]");
  }

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixFIRECuda::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE_CUDA;
  mask |= FINAL_INTEGRATE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

// fix_fire is a copy of fix_nve with modified velocity and time steps
// See Bitzek PRL 2006 
void FixFIRECuda::init()
{
  dtv = update->dt; // Initial time step is from LAMMPS
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv); // This would limit the step size if implemented

  cut = minsteps;
  mix = mix_in;

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"respa not supported");

  triggerneighsq = cuda->shared_data.atom.triggerneighsq;
  cuda->neighbor_decide_by_integrator = 1;
  Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);
  timesincelastdtupload=0;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

// Nearly identical to the initial_integrate of fix_nve_cuda.cpp
void FixFIRECuda::initial_integrate(int vflag)
{ 
// Evidently: If new, update this object's copy and also device's NVE object copy 
// of triggerneighsq (maximum sq dist before reneigh)
// Also update dtv, and dtf on the device
  if(triggerneighsq != cuda->shared_data.atom.triggerneighsq) {
    triggerneighsq = cuda->shared_data.atom.triggerneighsq;
    Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);
    timesincelastdtupload=0;
  }

  // The following variables would limit dx and prevent blow up.
  // Would need to implement nve/limit/cuda
  double dtfm,vsq,scale; 

  // update v and x of atoms in group

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  Cuda_FixNVECuda_InitialIntegrate(&cuda->shared_data, groupbit, nlocal);
}

/* ---------------------------------------------------------------------- */

// Similar to final_integrate of fix_nve_cuda.cpp
void FixFIRECuda::final_integrate()
{
  double dtfm,vsq,scale;

  // update v of atoms in group

  double *dev_v = static_cast<double*>(cuda->cu_v->dev_data());
  double *dev_f = static_cast<double*>(cuda->cu_f->dev_data());
  int *dev_mask = static_cast<int*>(cuda->cu_mask->dev_data());
  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  Cuda_FixNVECuda_FinalIntegrate(&cuda->shared_data, groupbit,nlocal);

  // Turn the global velocity a little bit  more along the global force...

  double vf1, vg_dot_vg1, Fg_dot_Fg1;

  // Later: This analyis could be done on the device, minimizing the number of to-host transfers
  // However, if certain fixes need the updated time step, more care will be needed.

  // Step F1 in Bitzek PRL 2006
  CU_CHECK( fix_fire_cuda_dot_products(groupbit, nlocal, nmax, dev_mask, dev_v,
				       dev_f, vf1, vg_dot_vg1, Fg_dot_Fg1) );

  double vf, vg_dot_vg, Fg_dot_Fg;
  MPI_Allreduce(&vf1, &vf, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&vg_dot_vg1, &vg_dot_vg, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&Fg_dot_Fg1, &Fg_dot_Fg, 1, MPI_DOUBLE, MPI_SUM, world);

  //  printf("%i, vf1 = %e, vf = %e\n", comm->me, vf1, vf);
  //  printf("%i, vg_dot_vg1 = %f, vg = %f\n", comm->me, vg_dot_vg1, vg_dot_vg1);
  //  printf("%i, Fg = %f, Fg = %f\n", comm->me, Fg_dot_Fg1, Fg_dot_Fg1);

  double help = 0.0;
  if (Fg_dot_Fg > 0.0) {
    help = mix*sqrt(vg_dot_vg/Fg_dot_Fg);
  }

  // First parameter should be (1-alpha), second should be alpha*|v|/|f|
  CU_CHECK( fix_fire_cuda_mix(1.0-mix, help, groupbit, nlocal, nmax, dev_mask,
			      dev_v, dev_f) );

  // Cut the velocities if the total power done by forces is negative
  // We prefer to have vf=0.0 (numerically) treated as the second case
  if (vf < 0.0) {
    CU_CHECK( fix_fire_cuda_mix(0.0, 0.0, groupbit, nlocal, nmax, dev_mask,
				dev_v, dev_f) );
    cut = minsteps;
    dtv = dtv*decfac;
    mix = mix_in;
  }
  else {
    // mixing is important only right after cut...
    if (cut < 0) {
      dtv = MIN(dtv*incfac, max_dt);
      mix = mix*mixdec;
    }
    else {
      cut = cut - 1;
    }
  }
  
  timesincelastdtupload++;
  if (dtv != update->dt) { // Must be careful if going to use this, because this fix's dt will be out-of-sync with the "true" (=everywhere else's) dt. && timesincelastdtupload > delayuploadfrequency) {
    timesincelastdtupload=0;
    update->dt = dtv;
    // For CUDA fixes, this updates the time step on the device, too.
    for (int i = 0; i < modify->nfix; i++) modify->fix[i]->reset_dt();
  }
}

/* ---------------------------------------------------------------------- */

void FixFIRECuda::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  Cuda_FixNVECuda_Init(&cuda->shared_data,dtv,dtf);
  // Every fix that needs the time step must update it on the GPU with reset_dt()
  
}
