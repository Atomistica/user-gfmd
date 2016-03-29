/*
   This software is distributed under the GNU General Public License.
Developing from fix_viscous_cuda.cu Tristan Sharp
------------------------------------------------------------------------- */

#include <stdio.h>
#define MY_PREFIX fix_addforce_dynamic_cuda
#include "cuda_shared.h"
#include "cuda_common.h"

#include "crm_cuda_utils.cu"

#include "fix_addforce_dynamic_cuda_cu.h"
#include "fix_addforce_dynamic_cuda_kernel.cu"

// Since we are not doing the force reduce in shared memory (yet) we dont need to
// Increase the size of sdata->buffer to hold reduces.
//void Cuda_FixAddForce_DynamicCuda_UpdateBuffer(cuda_shared_data* sdata)
//{
//  int3 layout = getgrid(sdata->atom.nlocal, 4 * sizeof(F_FLOAT));
//  dim3 threads(layout.z, 1, 1);
//  dim3 grid(layout.x, layout.y, 1);
//  int size = (unsigned)(layout.z * layout.y * layout.x) * 4 * sizeof(F_FLOAT);
//
//  if(sdata->buffersize < size) {
//    MYDBG(printf("Cuda_FixAddForce_DyanmicCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
//    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
//    sdata->buffer = CudaWrapper_AllocCudaData(size);
//    sdata->buffersize = size;
//    sdata->buffer_new++;
//    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
//  }
//
//  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
//}

void Cuda_FixAddForce_DynamicCuda_UpdateNmax(cuda_shared_data* sdata)
{
  // Update fix_afd pointers to point to the authoritative copy
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)    , & sdata->atom.type .dev_data, sizeof(int*));
}

void Cuda_FixAddForce_DynamicCuda_Init(cuda_shared_data* sdata)
{
  Cuda_FixAddForce_DynamicCuda_UpdateNmax(sdata);

}


void Cuda_FixAddForce_DynamicCuda_PostForce(cuda_shared_data* sdata, int groupbit,
		    double reflxcenter_sigma, double reflycenter_sigma,
	            double radspersigma_x, double radspersigma_y,
		    double acontact_sigma, 
		    double border_dist,
		    double r_overEstar)
{
  if(sdata->atom.update_nmax)
    Cuda_FixAddForce_DynamicCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  // Unlike fix_addforce_cuda, dynamic version does not compute scalar nor compute_vector so does not need buffer(?)
  //if(sdata->buffer_new)
  //  Cuda_FixAddForceCuda_UpdateBuffer(sdata);

  //int3 layout = getgrid(sdata->atom.nlocal, 4 * sizeof(F_FLOAT));
  int3 layout = getgrid(sdata->atom.nlocal, 0);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  //            bytes dynamically allocated per block of shared (used by vars declared as extern array)
  //Cuda_FixAddForce_DynamicCuda_PostForce_Kernel <<< grid, threads, threads.x* 3* sizeof(F_FLOAT)>>> (groupbit, 
  Cuda_FixAddForce_DynamicCuda_PostForce_Kernel <<< grid, threads>>> (groupbit, 
		    reflxcenter_sigma, reflycenter_sigma,
	            radspersigma_x, radspersigma_y,
		    acontact_sigma, 
		    border_dist,
		    r_overEstar);

  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Cuda_FixAddForce_DynamicCuda_PostForce: Kernel execution failed");

  //int oldgrid = grid.x;
  //grid.x = 4;
  //threads.x = 512;
  //reduce_foriginal <<< grid, threads, threads.x* sizeof(F_FLOAT)>>> (oldgrid, aforiginal);
  //cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Cuda_FixAddForce_DynamicCuda_PostForce: Kernel execution failed");
  // TAS Compute vector will here require that Reduce call maybe

}
