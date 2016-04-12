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
#include <stdio.h>

#include "cuda_common.h"

#include "fix_gfmd_cuda_cu.h"

#include "fix_gfmd_cuda_kernel.cu"

#define THREADS_PER_BLOCK 256

#define CU_CHECK(x)  { cudaError err = x; if (err != cudaSuccess) return err; }


cudaError fix_gfmd_cuda_list_to_grid(double xprd, double yprd, 
				     double xprd_half, double yprd_half,
				     int groupbit, int nlocal, int nall,
				     int nmax, int *dev_mask, X_FLOAT *dev_x,
				     int *dev_gid, X_FLOAT *dev_xeq,
				     int dxshift, int dyshift,
				     int nx, int ny, int xlo_loc, int ylo_loc,
				     int nx_loc, int ny_loc, int nxy_loc,
				     double *dev_u_xy, int &natoms_cur,
				     int &ngfmd_loc)
{
  int num_blocks = nall / THREADS_PER_BLOCK;
  //printf("list to grid: displacements fromdev %p todev %p \n",dev_x,dev_u_xy);//TAS

  if (num_blocks*THREADS_PER_BLOCK < nall)
    num_blocks++;

  int natoms_cur_buf[num_blocks], ngfmd_loc_buf[num_blocks];

  int *dev_natoms_cur_buf, *dev_ngfmd_loc_buf;
  CU_CHECK( cudaMalloc(&dev_natoms_cur_buf, num_blocks*sizeof(int)) );//printf("#CUDA Allocgfn %ubytes at%p\n",num_blocks*sizeof(int), dev_natoms_cur_buf);
  CU_CHECK( cudaMalloc(&dev_ngfmd_loc_buf, num_blocks*sizeof(int)) );//printf("#CUDA Allocgfn %ubytes at%p\n",num_blocks*sizeof(int), dev_ngfmd_loc_buf);

  fix_gfmd_cuda_list_to_grid_kernel
      <<<num_blocks, THREADS_PER_BLOCK, 2*THREADS_PER_BLOCK*sizeof(int)>>>
    (xprd, yprd, xprd_half, yprd_half, groupbit, nlocal, nall, nmax, dev_mask,
     dev_x, dev_gid, dev_xeq, dxshift, dyshift, nx, ny, xlo_loc, ylo_loc,
     nx_loc, ny_loc, nxy_loc, dev_u_xy, dev_natoms_cur_buf,
     dev_ngfmd_loc_buf);
  CU_CHECK( cudaPeekAtLastError() );

  CU_CHECK( cudaMemcpy(natoms_cur_buf, dev_natoms_cur_buf,
		       num_blocks*sizeof(int), cudaMemcpyDeviceToHost) );//printf("#CUDA CPgf to host %u bytes\n",num_blocks*sizeof(int));//TAS
  CU_CHECK( cudaMemcpy(ngfmd_loc_buf, dev_ngfmd_loc_buf,
		       num_blocks*sizeof(int), cudaMemcpyDeviceToHost) );//printf("#CUDA CPgf to host %u bytes\n",num_blocks*sizeof(int));//TAS

  for (int j = 0; j < num_blocks; j++) {
    natoms_cur += natoms_cur_buf[j];
    ngfmd_loc += ngfmd_loc_buf[j];
  }

  CU_CHECK( cudaFree(dev_natoms_cur_buf) );//printf("#CUDA free %p \n",dev_natoms_cur_buf);//TAS
  CU_CHECK( cudaFree(dev_ngfmd_loc_buf) );//printf("#CUDA free %p \n",dev_ngfmd_loc_buf);//TAS

  return cudaSuccess;
}



cudaError fix_gfmd_cuda_grid_to_list(int groupbit, int nlocal, int nall,
				     int nmax, int *dev_mask, F_FLOAT *dev_data,
				     int *dev_gid, int nx_loc, int ny_loc,
				     int nxy_loc, int xlo_loc, int ylo_loc,
				     int xhi_loc, int yhi_loc, double *dev_grid,
				     bool *dev_force_is_known, int &natoms_cur,
				     double sum_loc[3])
{
  int num_blocks = nall / THREADS_PER_BLOCK;

  //printf("on device grid to list: displacements fromdev %p todev %p \n",dev_grid,dev_data); // TAS

  if (num_blocks*THREADS_PER_BLOCK < nall)
    num_blocks++;

  int natoms_cur_buf[num_blocks];
  double sum_loc_x_buf[num_blocks], sum_loc_y_buf[num_blocks];
  double sum_loc_z_buf[num_blocks];

  int *dev_natoms_cur_buf;
  double *dev_sum_loc_x_buf, *dev_sum_loc_y_buf, *dev_sum_loc_z_buf;
  cudaMalloc(&dev_natoms_cur_buf, num_blocks*sizeof(int));//printf("#CUDA Allocgfngtl %ubytes at%p\n",num_blocks*sizeof(int), dev_natoms_cur_buf);
  cudaMalloc(&dev_sum_loc_x_buf, num_blocks*sizeof(double));//printf("#CUDA Allocgfngtl %ubytes at%p\n",num_blocks*sizeof(double), dev_sum_loc_x_buf);
  cudaMalloc(&dev_sum_loc_y_buf, num_blocks*sizeof(double));//printf("#CUDA Allocgfngtl %ubytes at%p\n",num_blocks*sizeof(double), dev_sum_loc_y_buf);
  cudaMalloc(&dev_sum_loc_z_buf, num_blocks*sizeof(double));//printf("#CUDA Allocgfngtl %ubytes at%p\n",num_blocks*sizeof(double), dev_sum_loc_z_buf);

  // going to add from grid to data: dev_grid to dev_data 
  fix_gfmd_cuda_grid_to_list_kernel
      <<<num_blocks, THREADS_PER_BLOCK, 4*THREADS_PER_BLOCK*sizeof(double)>>>
    (groupbit, nlocal, nall, nmax, dev_mask, dev_data, dev_gid, nx_loc, ny_loc,
     nxy_loc, xlo_loc, ylo_loc, xhi_loc, yhi_loc, dev_grid, dev_force_is_known,
     dev_natoms_cur_buf, dev_sum_loc_x_buf, dev_sum_loc_y_buf,
     dev_sum_loc_z_buf);
  cudaError err = cudaPeekAtLastError();
  if (err != cudaSuccess)
    return err;

  err = cudaMemcpy(natoms_cur_buf, dev_natoms_cur_buf,
		   num_blocks*sizeof(int), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    return err;
  }
  err = cudaMemcpy(sum_loc_x_buf, dev_sum_loc_x_buf, num_blocks*sizeof(double),
		   cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    return err;
  }
  err = cudaMemcpy(sum_loc_y_buf, dev_sum_loc_y_buf, num_blocks*sizeof(double),
		   cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    return err;
  }
  err = cudaMemcpy(sum_loc_z_buf, dev_sum_loc_z_buf, num_blocks*sizeof(double),
		   cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    return err;
  }

  natoms_cur = 0;
  for (int j = 0; j < num_blocks; j++) {
    natoms_cur += natoms_cur_buf[j];
    sum_loc[0] += sum_loc_x_buf[j];
    sum_loc[1] += sum_loc_y_buf[j];
    sum_loc[2] += sum_loc_z_buf[j];
  }

  cudaFree(dev_natoms_cur_buf);//printf("free%p\n",dev_natoms_cur_buf);//TAS
  cudaFree(dev_sum_loc_x_buf);//printf("free%p\n",dev_sum_loc_x_buf);//TAS
  cudaFree(dev_sum_loc_y_buf);//printf("free%p\n",dev_sum_loc_y_buf);//TAS
  cudaFree(dev_sum_loc_z_buf);//printf("free%p\n",dev_sum_loc_z_buf);//TAS

  return cudaSuccess;
}



cudaError fix_gfmd_cuda_xcm(int groupbit, double xprd,
			    double yprd, double zprd, int nlocal, int nmax,
			    int *mask, int *type, double *mass, double *rmass,
			    double *x, tagint *image,
			    double cmone[3])
{
  int num_blocks = nlocal / THREADS_PER_BLOCK;

  if (num_blocks*THREADS_PER_BLOCK < nlocal)
    num_blocks++;

  double cx_out[num_blocks], cy_out[num_blocks], cz_out[num_blocks];
  double *dev_cx_out, *dev_cy_out, *dev_cz_out;

  CU_CHECK( cudaMalloc(&dev_cx_out, num_blocks*sizeof(double)) ); //printf("Allocations, copies, frees without printfs here\n");//TAS
  CU_CHECK( cudaMalloc(&dev_cy_out, num_blocks*sizeof(double)) );
  CU_CHECK( cudaMalloc(&dev_cz_out, num_blocks*sizeof(double)) ); //printf("pointers to locations %p %p %p\n",dev_cx_out,dev_cy_out,dev_cz_out);

  if (rmass) {
    fix_gfmd_cuda_xcm_kernel1
      <<<num_blocks, THREADS_PER_BLOCK, 3*THREADS_PER_BLOCK*sizeof(double)>>>
      (groupbit, xprd, yprd, zprd, nlocal, nmax, mask, rmass, x, image,
       dev_cx_out, dev_cy_out, dev_cz_out);
    if (cudaPeekAtLastError() != cudaSuccess)
      return cudaPeekAtLastError();
  } else {
    fix_gfmd_cuda_xcm_kernel2
      <<<num_blocks, THREADS_PER_BLOCK, 3*THREADS_PER_BLOCK*sizeof(double)>>>
      (groupbit, xprd, yprd, zprd, nlocal, nmax, mask, type, mass, x, image,
       dev_cx_out, dev_cy_out, dev_cz_out);
    if (cudaPeekAtLastError() != cudaSuccess)
      return cudaPeekAtLastError();
  }

  CU_CHECK( cudaMemcpy(cx_out, dev_cx_out, num_blocks*sizeof(double),
		       cudaMemcpyDeviceToHost) );//printf("Copy to host x3\n");//TAS
  CU_CHECK( cudaMemcpy(cy_out, dev_cy_out, num_blocks*sizeof(double),
		       cudaMemcpyDeviceToHost) );
  CU_CHECK( cudaMemcpy(cz_out, dev_cz_out, num_blocks*sizeof(double),
		       cudaMemcpyDeviceToHost) );

  CU_CHECK( cudaFree(dev_cx_out) );//printf("freeing locations %p %p %p\n",dev_cx_out,dev_cy_out,dev_cz_out);
  CU_CHECK( cudaFree(dev_cy_out) );
  CU_CHECK( cudaFree(dev_cz_out) );

  cmone[0] = cmone[1] = cmone[2] = 0.0;
  for (int i = 0; i < num_blocks; i++) {
    cmone[0] += cx_out[i];
    cmone[1] += cy_out[i];
    cmone[2] += cz_out[i];
  }

  return cudaSuccess;
}


cudaError CudaData_printmyaddress(void* device_ptr) {  printmyaddress<<<1,1>>>(device_ptr); return cudaSuccess;} // TAS
//void CudaData_printmyaddre(float* device_ptr) {  printmyaddre<<<1,1>>>(device_ptr); } // TAS
//void CudaData_printmyaddre(double* device_ptr) {  printmyaddre<<<1,1>>>(device_ptr); } // TAS

