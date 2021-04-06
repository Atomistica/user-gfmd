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
