/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(gfmd/cuda,FixGFMDCuda)

#else

#ifndef FIX_GFMD_CUDA_H
#define FIX_GFMD_CUDA_H

#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "cuda_data.h"

#include "fix_gfmd.h"
	 
namespace LAMMPS_NS {

class FixGFMDCuda : public FixGFMD {

 public:
  FixGFMDCuda(class LAMMPS *, int, char **);
  virtual ~FixGFMDCuda();

  virtual int  setmask();
  virtual void init();
		virtual void setup(int);
  virtual void pre_force(int);
  virtual void post_force(int);

 protected:
  dev_array dev_u_xy;
  cCudaData<double, X_FLOAT, xy>* cu_u_xy;
  dev_array dev_f_xy;
  cCudaData<double, X_FLOAT, xy>* cu_f_xy;

  dev_array dev_force_is_known;
  cCudaData<bool, bool, x>* cu_force_is_known;

  class Cuda *cuda;              // pointer to cuda device data structures
    
  void grid_to_list_cuda();

  void xcm(int, double, double *);

};

}

#endif
#endif
