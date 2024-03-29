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

#include "cuda_shared.h"     // nve has it
#include "cuda_common.h"
#include "crm_cuda_utils.cu" // nve has it

#include "fix_fire_cuda_cu.h"
#include "fix_fire_cuda_kernel.cu"

#define THREADS_PER_BLOCK 256

#define CU_CHECK(x)  { cudaError err = x; if (err != cudaSuccess) return err; }


/*!
 * Launcher for array reduction
 */
template<typename T>
cudaError reduce(int n, T *arr, T &res)
{
  int ncur = n;

  int num_blocks = ncur / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < ncur)
    num_blocks++;

  T *buf1, *buf2;
  // CudaWrapper_AllocCudaData(unsigned nbytes); // Here we allocate memory. report to developer in printf // TAS
  CU_CHECK( cudaMalloc(&buf1, num_blocks*sizeof(T)) ); //printf("#CUDA Allocbuf1 %ubytes at%p\n",num_blocks*sizeof(T), buf1);
  CU_CHECK( cudaMalloc(&buf2, num_blocks*sizeof(T)) );//printf("#CUDA Allocbuf2 %ubytes at%p\n",num_blocks*sizeof(T), buf2);

  /*
   * First reduction step. Now have num_blocks entries in buf2
   */
  reduce_kernel<T, THREADS_PER_BLOCK>
    <<<num_blocks, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(T)>>>
    (ncur, arr, buf2);
  CU_CHECK( cudaPeekAtLastError() );

  /*
   * We need to reduce these num_blocks entries
   */
  ncur = num_blocks;
  while (ncur > 1) {
    /*
     * Flip buf1 and buf2. buf2 becomes input buffer
     */
    T *tmp = buf1;
    buf1 = buf2;
    buf2 = tmp;

    num_blocks = ncur / THREADS_PER_BLOCK;
    if (num_blocks*THREADS_PER_BLOCK < ncur)
      num_blocks++;

    /*
     * Reduction step.
     */
    reduce_kernel<T, THREADS_PER_BLOCK>
      <<<num_blocks, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(T)>>>
      (ncur, buf1, buf2);
    CU_CHECK( cudaPeekAtLastError() );

    ncur = num_blocks;
  }
  // Here we copy to host.  Does it get reported? TAS
  CU_CHECK( cudaMemcpy(&res, buf2, sizeof(T), cudaMemcpyDeviceToHost) );//printf("#CUDA CP to host %u bytes\n",sizeof(T));//TAS

  CU_CHECK( cudaFree(buf1) );//printf("#CUDA free buf1%p \n",buf1);//TAS
  CU_CHECK( cudaFree(buf2) );//printf("#CUDA free buf2%p \n",buf2);//TAS

  return cudaSuccess;
}

/*!
 * Compute the dot products v.f, v.v and f.f
 */
cudaError fix_fire_cuda_dot_products(int groupbit, int nlocal, int nmax,
                                     int *dev_mask, double *dev_v,
                                     double *dev_f, double &vf,
                                     double &vg_dot_vg, double &Fg_dot_Fg)
{
  int num_blocks = nlocal / THREADS_PER_BLOCK;

  if (num_blocks*THREADS_PER_BLOCK < nlocal)
    num_blocks++;

  double *dev_vf, *dev_vg_dot_vg, *dev_Fg_dot_Fg;
  CU_CHECK( cudaMalloc(&dev_vf, num_blocks*sizeof(double)) );//printf("#CUDA alloc%ubytes at%p \n",num_blocks*sizeof(double), dev_vf);//TAS
  CU_CHECK( cudaMalloc(&dev_vg_dot_vg, num_blocks*sizeof(double)) );//printf("#CUDA alloc%ubytes at%p \n",num_blocks*sizeof(double), dev_vg_dot_vg);//TAS
  CU_CHECK( cudaMalloc(&dev_Fg_dot_Fg, num_blocks*sizeof(double)) );//printf("#CUDA alloc%ubytes at%p \n",num_blocks*sizeof(double), dev_Fg_dot_Fg);//TAS

  fix_fire_cuda_dot_products<THREADS_PER_BLOCK>
      <<<num_blocks, THREADS_PER_BLOCK, 3*THREADS_PER_BLOCK*sizeof(double)>>>
    (groupbit, nlocal, nmax, dev_mask, dev_v, dev_f, dev_vf, dev_vg_dot_vg,
     dev_Fg_dot_Fg);
  CU_CHECK( cudaPeekAtLastError() );

  // Reduce each thread's values and copy to host pointer
  CU_CHECK( reduce(num_blocks, dev_vf, vf) );
  CU_CHECK( reduce(num_blocks, dev_vg_dot_vg, vg_dot_vg) );
  CU_CHECK( reduce(num_blocks, dev_Fg_dot_Fg, Fg_dot_Fg) );

  CU_CHECK( cudaFree(dev_vf) );//printf("#CUDA free %p \n",dev_vf);//TAS
  CU_CHECK( cudaFree(dev_vg_dot_vg) );//printf("#CUDA free %p \n",dev_vg_dot_vg);//TAS
  CU_CHECK( cudaFree(dev_Fg_dot_Fg) );//printf("#CUDA free %p \n",dev_Fg_dot_Fg);//TAS

  return cudaSuccess;
}

/*!
 * Mix each atom's force vector into its velocity vector
 */
cudaError fix_fire_cuda_mix(double a, double b, int groupbit, int nlocal,
                            int nmax, int *dev_mask, double *dev_v,
                            double *dev_f)
{
  int num_blocks = nlocal / THREADS_PER_BLOCK;

  if (num_blocks*THREADS_PER_BLOCK < nlocal)
    num_blocks++;

  fix_fire_cuda_mix_kernel<<<num_blocks, THREADS_PER_BLOCK>>>
    (a, b, groupbit, nlocal, nmax, dev_mask, dev_v, dev_f);
  CU_CHECK( cudaPeekAtLastError() );

  return cudaSuccess;
}
