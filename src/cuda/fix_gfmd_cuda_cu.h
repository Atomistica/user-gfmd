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
#ifndef FIX_GFMD_CUDA_CU_H
#define FIX_GFMD_CUDA_CU_H

#include <cuda_runtime.h>
//#include <mpi.h>

#include "../../src/lmptype.h"

using namespace LAMMPS_NS;

extern "C" {

cudaError fix_gfmd_cuda_list_to_grid(double xprd, double yprd, 
				     double xprd_half, double yprd_half,
				     int groupbit, int nlocal, int nall,
				     int nmax, int *mask, X_FLOAT *x, int *gid,
				     X_FLOAT *xeq, int dxshift, int dyshift,
				     int nx, int ny, int xlo_loc, int ylo_loc,
				     int nx_loc, int ny_loc, int nxy_loc,
				     double *u_xy, int &natoms_cur,
				     int &ngfmd_loc);

cudaError fix_gfmd_cuda_grid_to_list(int groupbit, int nlocal, int nall,
				     int nmax, int *dev_mask, F_FLOAT *dev_data,
				     int *dev_gid, int nx, int ny, int nxy,
				     int xlo_loc, int ylo_loc, int xhi_loc,
				     int yhi_loc, double *dev_grid,
				     bool *dev_force_is_known, int &natoms_cur,
				     double sum_loc[3]);

cudaError fix_gfmd_cuda_xcm(int groupbit, double xprd,
			    double yprd, double zprd, int nlocal, int nmax,
			    int *mask, int *type, double *mass, double *rmass,
			    double *x, tagint *image,
			    double cmone[3]);

cudaError CudaData_printmyaddress(void* device_ptr);// TAS
//void CudaData_printmyaddre(float* device_ptr);// TAS
//void CudaData_printmyaddre(double* device_ptr);

}

#endif
