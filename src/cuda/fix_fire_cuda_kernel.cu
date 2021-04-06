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
//#include <stdio.h>

//#include "cuda_common.h"

//#include "fix_fire_cuda_cu.h"

#define THREADS_PER_BLOCK 256

template<typename T, unsigned int blockDim_x>
__device__ void _warp_reduce(volatile T *sbuf, int tid)
{
  if (blockDim_x >= 64)  sbuf[tid] += sbuf[tid + 32];
  if (blockDim_x >= 32)  sbuf[tid] += sbuf[tid + 16];
  if (blockDim_x >= 16)  sbuf[tid] += sbuf[tid +  8];
  if (blockDim_x >= 8)   sbuf[tid] += sbuf[tid +  4];
  if (blockDim_x >= 4)   sbuf[tid] += sbuf[tid +  2];
  if (blockDim_x >= 2)   sbuf[tid] += sbuf[tid +  1];
}


/*!
 * Reduction of three scalar variables of the same type
 */
template<typename T, unsigned int blockDim_x>
__device__ void _reduce3(T x1, T *x1_out, T x2, T *x2_out, T x3, T *x3_out)
{
  extern __shared__ T sbuf[];

  sbuf[threadIdx.x] = x1;						
  sbuf[threadIdx.x+blockDim_x] = x2;					
  sbuf[threadIdx.x+2*blockDim_x] = x3;					
  __syncthreads();							
  
  for (int s = blockDim_x>>1; s > 32; s >>= 1) {			
    if (threadIdx.x < s) {						
      sbuf[threadIdx.x] += sbuf[threadIdx.x+s];				
      sbuf[threadIdx.x+blockDim_x] += sbuf[threadIdx.x+blockDim_x+s];	
      sbuf[threadIdx.x+2*blockDim_x] += sbuf[threadIdx.x+2*blockDim_x+s]; 
    }									
    __syncthreads();							
  }									

  if (threadIdx.x < 32) {						
    _warp_reduce<T, blockDim_x>(sbuf, threadIdx.x);
    _warp_reduce<T, blockDim_x>(sbuf, threadIdx.x+blockDim_x);
    _warp_reduce<T, blockDim_x>(sbuf, threadIdx.x+2*blockDim_x);
  }
  
  if (threadIdx.x == 0) {						
    x1_out[blockIdx.x] = sbuf[0];					
    x2_out[blockIdx.x] = sbuf[blockDim_x];				
    x3_out[blockIdx.x] = sbuf[2*blockDim_x];				
  }
}


/*!
 * Kernel for array reduction
 */
template<typename T, unsigned int blockDim_x>
__global__ void reduce_kernel(int n, T *idata, T *odata)
{
  extern __shared__ T sdata[];

  int i = blockIdx.x*blockDim_x + threadIdx.x;
  int grid_size = blockDim_x*gridDim.x;

  sdata[threadIdx.x] = 0.0;
  while (i < n) {
    sdata[threadIdx.x] += idata[i];
    i += grid_size;
  }
  __syncthreads();

  for (int s = blockDim_x>>1; s > 32; s >>= 1) {			
    if (threadIdx.x < s) {						
      sdata[threadIdx.x] += sdata[threadIdx.x+s];				
    }									
    __syncthreads();							
  }									

  if (threadIdx.x < 32) {
    _warp_reduce<T, blockDim_x>(sdata, threadIdx.x);
  }

  if (threadIdx.x == 0) {
    odata[blockIdx.x] = sdata[0];
  }
}



/*!
 * Compute the dot products v.f, v.v and f.f
 */
template<unsigned int blockDim_x>
__global__ void fix_fire_cuda_dot_products(int groupbit, int nlocal, int nmax,
					   int *mask, double *v, double *f,
					   double *vf, double *vg_dot_vg,
					   double *Fg_dot_Fg)
{
  int i = blockIdx.x*blockDim_x + threadIdx.x;

  /*
   * Compute dot products
   */

  double vf1, vg_dot_vg1, Fg_dot_Fg1;
  vf1 = vg_dot_vg1 = Fg_dot_Fg1 = 0.0;

  if (i < nlocal) {
    if (mask[i] & groupbit) {
      vf1        += v[i]*f[i] + v[i+nmax]*f[i+nmax] + v[i+2*nmax]*f[i+2*nmax];
      vg_dot_vg1 += v[i]*v[i] + v[i+nmax]*v[i+nmax] + v[i+2*nmax]*v[i+2*nmax];
      Fg_dot_Fg1 += f[i]*f[i] + f[i+nmax]*f[i+nmax] + f[i+2*nmax]*f[i+2*nmax];
    }
  }

  /*
   * Reduce counters
   */
  
  _reduce3<double, blockDim_x>
    (vf1, vf, vg_dot_vg1, vg_dot_vg, Fg_dot_Fg1, Fg_dot_Fg);
}


/*!
 * Mix velocities and forces
 */
__global__ void fix_fire_cuda_mix_kernel(double a, double b, int groupbit,
					 int nlocal, int nmax, int *mask,
					 double *v, double *f) // May be able to use atom->v later
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;

  /*
   * Compute dot products
   */

  if (i < nlocal) {
    if (mask[i] & groupbit) {
      v[i       ] = a*v[i       ] + b*f[i       ];
      v[i+nmax  ] = a*v[i+  nmax] + b*f[i+  nmax];
      v[i+2*nmax] = a*v[i+2*nmax] + b*f[i+2*nmax];
    }
  }
}


