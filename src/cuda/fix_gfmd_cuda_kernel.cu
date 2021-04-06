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
#include <stdio.h>

#include "cuda_common.h"

#include "../../src/gfmd_grid.h"

// #include "fix_gfmd_cuda_kernel.h"

/*
 * This should give the same answer as the Fortran MODULO intrinsic.
 * I.e. the answer of modulo(x, y) will always be within [0..y-1]
 */
#define modulo(x, y)  ( (x) < 0 ? (x)-(((x)-(y)+1)/(y))*(y) : (x)-((x)/(y))*(y) )
#define fmodulo(x, y)  ( (x)-floor((x)/(y))*(y) )
// Parentheses for order of operations

#define THREADS_PER_BLOCK 256

#define CU_CHECK(x)  { cudaError err = x; if (err != cudaSuccess) return err; }


template<typename T>
__device__ void warp_reduce(volatile T *sbuf, int tid)
{
  sbuf[tid] += sbuf[tid + 32];
  sbuf[tid] += sbuf[tid + 16];
  sbuf[tid] += sbuf[tid +  8];
  sbuf[tid] += sbuf[tid +  4];
  sbuf[tid] += sbuf[tid +  2];
  sbuf[tid] += sbuf[tid +  1];
}


template<typename T>
__device__ void reduce3(T x1, T *x1_out, T x2, T *x2_out, T x3, T *x3_out)
{
  extern __shared__ T sbuf[];

  sbuf[threadIdx.x] = x1;						
  sbuf[threadIdx.x+blockDim.x] = x2;					
  sbuf[threadIdx.x+2*blockDim.x] = x3;					
  __syncthreads();			     
  
  for (int s = blockDim.x>>1; s > 32; s >>= 1) {			
    if (threadIdx.x < s) {						
      sbuf[threadIdx.x] += sbuf[threadIdx.x+s];				
      sbuf[threadIdx.x+blockDim.x] += sbuf[threadIdx.x+blockDim.x+s];	
      sbuf[threadIdx.x+2*blockDim.x] += sbuf[threadIdx.x+2*blockDim.x+s]; 
    }									
    __syncthreads();							
  }									

  if (threadIdx.x < 32) {						
    warp_reduce(sbuf, threadIdx.x);					
    warp_reduce(sbuf, threadIdx.x+blockDim.x);				
    warp_reduce(sbuf, threadIdx.x+2*blockDim.x);
  }
  
  if (threadIdx.x == 0) {						
    x1_out[blockIdx.x] = sbuf[0];					
    x2_out[blockIdx.x] = sbuf[blockDim.x];				
    x3_out[blockIdx.x] = sbuf[2*blockDim.x];				
  }
}


__global__ void fix_gfmd_cuda_list_to_grid_kernel(double xprd, double yprd,
						  double xprd_half, 
						  double yprd_half,
						  int groupbit,
						  int nlocal, int nall,
						  int nmax, int *mask,
						  X_FLOAT *x, int *gid,
						  X_FLOAT *xeq,
						  int dxshift, int dyshift,
						  int nx, int ny, int xlo_loc,
						  int ylo_loc, int nx_loc,
						  int ny_loc, int nxy_loc,
						  double *u_xy,
						  int *natoms_cur_out,
						  int *ngfmd_loc_out)
{
  extern __shared__ int sdata[];

  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int natoms_cur = 0, ngfmd_loc = 0;

  // FIXME! i < nall does not work at the first time step. The ghost atoms
  // have wrong xeq/gid?
  //  if (i < nall) {
  if (i < nlocal) {
    if (mask[i] & groupbit) {
      int j = gid[i];

      int ix = IX_FROM_POW2_IDX(j) - dxshift;
      int iy = IY_FROM_POW2_IDX(j) - dyshift;
      int iu = IU_FROM_POW2_IDX(j);

      if (dxshift != 0 || dyshift != 0) {
	// wrap to box
	ix = modulo(ix, nx);
	iy = modulo(iy, ny);

	// store wrapped grid id for unpacking/communication
	gid[i] = POW2_IDX(ix, iy, iu);
      }

      ix -= xlo_loc;
      iy -= ylo_loc;

      if (ix >= 0 && ix < nx_loc && iy >= 0 && iy < ny_loc) {

	// y is the fast index
	int iloc = ix*ny_loc + iy;
	int idof = 3*iu;

	double ux = x[i] - xeq[i];
	double uy = x[i+nmax] - xeq[i+nmax];
	double uz = x[i+2*nmax] - xeq[i+2*nmax];

	// wrap to box
	ux = fmodulo(ux+xprd_half, xprd)-xprd_half;
	uy = fmodulo(uy+yprd_half, yprd)-yprd_half;

	u_xy[ idof   *nxy_loc+iloc] = ux;
	u_xy[(idof+1)*nxy_loc+iloc] = uy;
	u_xy[(idof+2)*nxy_loc+iloc] = uz;

	natoms_cur++;
      }

      if (i < nlocal) {
	ngfmd_loc++;
      }
    }
  }

  /*
   * Reduce sanity check counters natoms_cur and ngfmd_loc
   */
  sdata[threadIdx.x] = natoms_cur;
  sdata[threadIdx.x+blockDim.x] = ngfmd_loc;
  __syncthreads();

  for (int s = blockDim.x>>1; s > 32; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x+blockDim.x] += sdata[threadIdx.x+blockDim.x + s];
    }
    __syncthreads();
  }

  if (threadIdx.x < 32) {
    warp_reduce(sdata, threadIdx.x);
    warp_reduce(sdata, blockDim.x+threadIdx.x);
  }

  if (threadIdx.x == 0) {
    natoms_cur_out[blockIdx.x] = sdata[0];
    ngfmd_loc_out[blockIdx.x] = sdata[blockDim.x];
  }
}


__global__ void fix_gfmd_cuda_grid_to_list_kernel(int groupbit,
						  int nlocal, int nall,
						  int nmax, int *mask,
						  F_FLOAT *data, int *gid,
						  int nx_loc, int ny_loc,
						  int nxy_loc, int xlo_loc,
						  int ylo_loc, int xhi_loc,
						  int yhi_loc, double *grid,
						  bool *force_is_known,
						  int *natoms_cur_out,
						  double *sum_loc_x_out,
						  double *sum_loc_y_out,
						  double *sum_loc_z_out)
{
  /*
   * Both __shared__ arrays point to the same memory location.
   * FIXME!!! This assumes that int has a smaller size than double.
   */
  extern __shared__ int sidata[];
  extern __shared__ double sddata[];

  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int natoms_cur = 0;
  double sum_loc_x = 0.0, sum_loc_y = 0.0, sum_loc_z = 0.0;

  // FIXME! i < nall does not work at the first time step. The ghost atoms
  // have wrong xeq/gid?
  //  if (i < nall) {
  if (i < nlocal) {
    if (mask[i] & groupbit) {
      int j = gid[i];

      int ix = IX_FROM_POW2_IDX(j);
      int iy = IY_FROM_POW2_IDX(j);

      // has this atom been computed locally?
      if (ix >= xlo_loc && ix <= xhi_loc && iy >= ylo_loc && iy <= yhi_loc) {

	ix -= xlo_loc;
	iy -= ylo_loc;

        int iu = IU_FROM_POW2_IDX(j);

        // y is the fast index
        int iloc = ix*ny_loc + iy;
        int idof = 3*iu;

	// cuda version of gfmd adds gfmd forces directly to the lammps force list
        data[i       ]  += grid[ idof   *nxy_loc+iloc];
        data[i+  nmax]  += grid[(idof+1)*nxy_loc+iloc];
	data[i+2*nmax]  += grid[(idof+2)*nxy_loc+iloc];
#ifdef GFMD_DEBUG
	if (i==0) printf("GFMD CUDA: Atom 0 z-force from gfmd:%f, total:%f \n",
                         grid[(idof+2)*nxy_loc+iloc],data[i+2*nmax]);
#endif
        if (i < nlocal) {
          sum_loc_x += data[i       ];
          sum_loc_y += data[i+  nmax];
          sum_loc_z += data[i+2*nmax];
          natoms_cur++;
        }

        force_is_known[i] = true;
      }
    }
  }

  /*
   * Reduce counters natoms_cur and sum_loc_*
   */

  sidata[threadIdx.x] = natoms_cur;
  sddata[threadIdx.x+blockDim.x] = sum_loc_x;
  sddata[threadIdx.x+2*blockDim.x] = sum_loc_y;
  sddata[threadIdx.x+3*blockDim.x] = sum_loc_z;
  __syncthreads();

  for (int s = blockDim.x>>1; s > 32; s >>= 1) {
    if (threadIdx.x < s) {
      sidata[threadIdx.x] += sidata[threadIdx.x + s];
      sddata[threadIdx.x+blockDim.x] += sddata[threadIdx.x+blockDim.x + s];
      sddata[threadIdx.x+2*blockDim.x] += sddata[threadIdx.x+2*blockDim.x + s];
      sddata[threadIdx.x+3*blockDim.x] += sddata[threadIdx.x+3*blockDim.x + s];
    }
    __syncthreads();
  }

  if (threadIdx.x < 32) {
    warp_reduce(sidata, threadIdx.x);
    warp_reduce(sddata, blockDim.x+threadIdx.x);
    warp_reduce(sddata, blockDim.x+2*threadIdx.x);
    warp_reduce(sddata, blockDim.x+3*threadIdx.x);
  }

  if (threadIdx.x == 0) {
    natoms_cur_out[blockIdx.x] = sidata[0];
    sum_loc_x_out[blockIdx.x] = sddata[blockDim.x];
    sum_loc_y_out[blockIdx.x] = sddata[2*blockDim.x];
    sum_loc_z_out[blockIdx.x] = sddata[3*blockDim.x];
  }
}


__global__ void fix_gfmd_cuda_xcm_kernel1(int groupbit, double xprd,
					  double yprd, double zprd, int nlocal,
					  int nmax, int *mask, double *rmass,
					  double *x, tagint *image,
					  double *cx_out, double *cy_out,
					  double *cz_out)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;

  double cx = 0.0, cy = 0.0, cz = 0.0;

  if (i < nlocal) {
    if (mask[i] & groupbit) {
      double massone = rmass[i];

      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      
      double a = x[i] + xbox*xprd;
      double b = x[nmax+i] + ybox*yprd;
      double c = x[2*nmax+i] + zbox*zprd;

      cx += a * massone;
      cy += b * massone;
      cz += c * massone;
    }
  }

  /*
   * Reduce cm counters
   */
  
  reduce3(cx, cx_out, cy, cy_out, cz, cz_out);
}


__global__ void fix_gfmd_cuda_xcm_kernel2(int groupbit, double xprd,
					  double yprd, double zprd, int nlocal,
					  int nmax, int *mask, int *type,
					  double *mass, double *x,
					  tagint *image,  double *cx_out,
					  double *cy_out, double *cz_out)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;

  double cx = 0.0, cy = 0.0, cz = 0.0;

  if (i < nlocal) {
    if (mask[i] & groupbit) {
      double massone = mass[type[i]];

      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      
      double a = x[i] + xbox*xprd;
      double b = x[nmax+i] + ybox*yprd;
      double c = x[2*nmax+i] + zbox*zprd;

      cx += a * massone;
      cy += b * massone;
      cz += c * massone;
    }
  }

  /*
   * Reduce cm counters
   */

  reduce3(cx, cx_out, cy, cy_out, cz, cz_out);
}



// DEBUG
__global__ void printmyaddress(void* ptr) {
  double* dptr = (double *)(ptr);
  float* fptr = (float *)(ptr);
           
  printf(" void ptr at device memory %p + 0 is dbl %f, or is flt %f)\n",ptr,dptr[0],fptr[0]);
  printf(" void ptr at device memory %p + 1 is dbl %f, or is flt %f)\n",ptr,dptr[1],fptr[1]);
  printf(" void ptr at device memory %p + 2 is dbl %f, or is flt %f)\n",ptr,dptr[2],fptr[2]);
  printf(" void ptr at device memory %p + 3 is dbl %f, or is flt %f)\n",ptr,dptr[3],fptr[3]);
  printf(" void ptr at device memory %p + 4 is dbl %f, or is flt %f)\n",ptr,dptr[4],fptr[4]);
  printf(" void ptr at device memory %p + 5 is dbl %f, or is flt %f)\n",ptr,dptr[5],fptr[5]);
  printf(" void ptr at device memory %p + 6 is dbl %f, or is flt %f)\n",ptr,dptr[6],fptr[6]);
  printf(" void ptr at device memory %p + 7 is dbl %f, or is flt %f)\n",ptr,dptr[7],fptr[7]);
  printf(" void ptr at device memory %p + 8 is dbl %f, or is flt %f)\n",ptr,dptr[8],fptr[8]);
  printf(" void ptr at device memory %p + 9 is dbl %f, or is flt %f)\n",ptr,dptr[9],fptr[9]);
  printf(" void ptr at device memory %p + 10 is dbl %f, or is flt %f)\n",ptr,dptr[10],fptr[10]);
  printf(" void ptr at device memory %p + 11 is dbl %f, or is flt %f)\n",ptr,dptr[11],fptr[11]);
  printf(" void ptr at device memory %p + 12 is dbl %f, or is flt %f)\n",ptr,dptr[12],fptr[12]);
  printf(" void ptr at device memory %p + 13 is dbl %f, or is flt %f)\n",ptr,dptr[13],fptr[13]);
  printf(" void ptr at device memory %p + 14 is dbl %f, or is flt %f)\n",ptr,dptr[14],fptr[14]);
  printf(" void ptr at device memory %p + 15 is dbl %f, or is flt %f)\n",ptr,dptr[15],fptr[15]);
  printf(" void ptr at device memory %p + 16 is dbl %f, or is flt %f)\n",ptr,dptr[16],fptr[16]);
  printf(" void ptr at device memory %p + 17 is dbl %f, or is flt %f)\n",ptr,dptr[17],fptr[17]);

  ptr = (void *)0x2300707000;
  dptr = (double *)ptr;
  fptr = (float *)ptr;
  printf("         special u: %p dbl %f or flt %f)\n",ptr,dptr[0],fptr[0]);
  ptr = (void *)0x2300707600;
  dptr = (double *)ptr;
  fptr = (float *)ptr;
  printf("         special f: %p dbl %f or flt %f)\n",ptr,dptr[0],fptr[0]);
}



