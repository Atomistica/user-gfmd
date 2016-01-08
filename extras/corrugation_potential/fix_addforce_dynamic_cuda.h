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

// TAS developing by referencing fix_viscous_cuda.h and fix_set_force_cuda.h

#ifdef FIX_CLASS

FixStyle(addforce/dynamic/cuda,FixAddForceDynamicCuda)

#else

#ifndef LMP_FIX_ADDFORCE_DYNAMIC_CUDA_H
#define LMP_FIX_ADDFORCE_DYNAMIC_CUDA_H

#include "fix_addforce_dynamic.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class FixAddForceDynamicCuda : public FixAddForceDynamic {
 public:
  FixAddForceDynamicCuda(class LAMMPS *, int, char **);
  ~FixAddForceDynamicCuda(); // Is this necessary to be called?
  int setmask();
  void setup(int); 
  //void min_setup(int);
  void post_force(int);
  //cCudaData<double, F_FLOAT, x>* cu_reflxcenter_sigma;
  //cCudaData<double, F_FLOAT, x>* cu_reflycenter_sigma;
  //cCudaData<double, F_FLOAT, x>* cu_radspersigma_x;
  //cCudaData<double, F_FLOAT, x>* cu_radspersigma_y;

 private:
  class Cuda *cuda;

  // We inherit what we need from fix_addforce_dynamic

  //double shiftx_nn_; // shift of hypothetical opposing lattice, units nn
  //double shifty_nn_; // shift of hypothetical opposing lattice, units nn
  //double acontact_nn_;  // effective contact radius, units of nn spacing
  //double nndist_;
  //double r_overEstar;
  
  //int force_flag;
  //ouble foriginal_all[4];
  //double foriginal[4];
 
  //  class Cuda *cuda;
  //  double cu_reflxcenter_sigma;
  //  double cu_reflycenter_sigma;
  //  double cu_radspersigma_x;
  //  double acontact_nn_;  // effective contact radius, units of nn spacing
  //  double nndist_;
  //  double r_overEstar;  double cu_radspersigma_y;
  //  
  //  int force_flag;
  //  double foriginal_all[4];
  //  double foriginal[4];
  };


}

#endif
#endif
