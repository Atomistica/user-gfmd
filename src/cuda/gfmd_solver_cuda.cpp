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
#include <string.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>

#include "gfmd_solver_cuda.h"

#include "pointers.h"

#include "comm.h"
#include "domain.h"
#include "linearalgebra.h"
#include "memory.h"
#include "mpi.h"
//#include "timer.h"

#include "gfmd_misc.h"
#include "gfmd_solver_cuda_cu.h"

using namespace LAMMPS_NS;

//#define AVOID_BATCH_TRANSFORM

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

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

/* ----------------------------------------------------------------------
 * CUDA helper stuff
 * --------------------------------------------------------------------*/

/*!
 * Check the CUDA memory type for this pointer
 */
inline cudaMemoryType cu_memory_type(void *data) {
  cudaPointerAttributes attr;
  cudaError_t err = cudaPointerGetAttributes(&attr, data);
  /* cudaPointerGetAttributes seems to generate an error for malloc'd ptrs */
  if (err == cudaErrorInvalidValue) {
    /* Clear error state */
    cudaGetLastError();
    return cudaMemoryTypeHost;
  }
  return attr.memoryType;
}

/*!
 * Check whether this matrix resides on the device
 */
inline bool cu_on_device(void *data) {
  return cu_memory_type(data) == cudaMemoryTypeDevice;
}

/*!
 * Check whether this matrix resides on the device
 */
inline bool cu_on_host(void *data) {
  return cu_memory_type(data) == cudaMemoryTypeHost;
}



/* ----------------------------------------------------------------------
 * Add some documentation
 * --------------------------------------------------------------------*/

GFMDSolverStaticCUDA::GFMDSolverStaticCUDA(LAMMPS *lmp) : GFMDSolver(lmp)
{
  strcpy(name, "static/cuda");

  linf_ = NULL;
  dev_linf_ = NULL;

  fft_plan_forward_ = NULL;
  fft_plan_inverse_ = NULL;

  dev_q_operator_ = NULL;
  dev_q_buffer1_ = NULL;
  dev_q_buffer2_ = NULL;

  dev_epot_buf_ = NULL;
}


GFMDSolverStaticCUDA::~GFMDSolverStaticCUDA()
{
  if (linf_)
    delete [] linf_;
  if (dev_linf_)
    cudaFree(dev_linf_);

  if (fft_plan_forward_)
    cufftDestroy(fft_plan_forward_);
  if (fft_plan_inverse_)
    cufftDestroy(fft_plan_inverse_);

  if (dev_q_operator_)
    cudaFree(dev_q_operator_);
  if (dev_q_buffer1_)
    cudaFree(dev_q_buffer1_);
  if (dev_q_buffer2_)
    cudaFree(dev_q_buffer2_);

  if (dev_epot_buf_)
    cudaFree(dev_epot_buf_);
}


void GFMDSolverStaticCUDA::init()
{
  size_t pitch;

  dev_epot_buf_ = NULL;

  nbuffer_half_loc_ = nx*(ny/2+1);

  if (nprocs != 1) {
    error->all(FLERR,"[GFMDSolverStaticCUDA::GFMDSolverStaticCUDA] CUDA FFT "
	       "for GFMD currently only works on a single processor");
  }

  /* Don't know if this is always true */
  if (sizeof(cufftDoubleReal) != sizeof(double)) {
    error->all(FLERR,"cufftDoubleReal and double size mismatch.");
  }
  if (sizeof(cufftDoubleComplex) != sizeof(double_complex)) {
    error->all(FLERR,"cufftDoubleComplex and double_complex size mismatch.");
  }

  /*
  if (cudaMallocPitch(&dev_q_operator_, &dev_q_operator_pitch_,
		      nbuffer_half_loc_*sizeof(cufftDoubleComplex), ndof_sq)
      != cudaSuccess) {
    error->all(FLERR,"cudaMalloc failed.");
  }
  */
  dev_q_operator_pitch_ = ndof_sq*sizeof(cufftDoubleComplex);
  if (dev_q_operator_ != NULL) cudaFree(dev_q_operator_);
  if (cudaMalloc(&dev_q_operator_, 
                 nbuffer_half_loc_*ndof_sq*sizeof(cufftDoubleComplex))
      != cudaSuccess) {
    error->all(FLERR,"cudaMalloc failed.");
  }

  if (dev_q_buffer1_ != NULL) cudaFree(dev_q_buffer1_);
  if (cudaMallocPitch(&dev_q_buffer1_, &dev_q_buffer_pitch_,
		      nbuffer_half_loc_*sizeof(cufftDoubleComplex), ndof)
      != cudaSuccess) {
    error->all(FLERR,"cudaMalloc failed.");
  }

  if (dev_q_buffer2_ != NULL) cudaFree(dev_q_buffer2_);
  if (cudaMallocPitch(&dev_q_buffer2_, &pitch,
		      nbuffer_half_loc_*sizeof(cufftDoubleComplex), ndof)
      != cudaSuccess) {
    error->all(FLERR,"cudaMalloc failed.");
  }

  if (pitch != dev_q_buffer_pitch_) {
    error->all(FLERR,"pitch != dev_q_buffer_pitch");
  }

  pitch = dev_q_operator_pitch_/sizeof(cufftDoubleComplex);
  if (pitch*sizeof(cufftDoubleComplex) != dev_q_operator_pitch_) {
    error->all(FLERR,"pitch*sizeof(cufftDoubleComplex) != "
	       "dev_q_operator_pitch");
  }
  dev_q_operator_pitch_ = pitch;

  pitch = dev_q_buffer_pitch_/sizeof(cufftDoubleComplex);
  if (pitch*sizeof(cufftDoubleComplex) != dev_q_buffer_pitch_) {
    error->all(FLERR,"pitch*sizeof(cufftDoubleComplex) != dev_q_buffer_pitch");
  }
  dev_q_buffer_pitch_ = pitch;

  if (fft_plan_forward_)
    cufftDestroy(fft_plan_forward_);
  if (fft_plan_inverse_)
    cufftDestroy(fft_plan_inverse_);

#ifdef AVOID_BATCH_TRANSFORM
  if (cufftPlan2d(&fft_plan_forward_, nx, ny, CUFFT_D2Z) != CUFFT_SUCCESS) {
    error->all(FLERR,"Failed to plan forward FFT.");
  }
  if (cufftPlan2d(&fft_plan_inverse_, nx, ny, CUFFT_Z2D) != CUFFT_SUCCESS) {
    error->all(FLERR,"Failed to plan inverse FFT.");
  }
#else
  int dims1[2] = { nx, ny };
  int dims2[2] = { nx, ny/2+1 };
  if (cufftPlanMany(&fft_plan_forward_, 2, dims1, dims1, 1, nxy_loc,
		    dims2, 1, dev_q_buffer_pitch_, CUFFT_D2Z, ndof)
      != CUFFT_SUCCESS) {
    error->one(FLERR,"Failed to plan forward FFT.");
  }
  if (cufftPlanMany(&fft_plan_inverse_, 2, dims1, dims2, 1, dev_q_buffer_pitch_,
		    dims1, 1, nxy_loc, CUFFT_Z2D, ndof)
      != CUFFT_SUCCESS) {
    error->one(FLERR,"Failed to plan forward FFT.");
  }
#endif

  if (cufftSetCompatibilityMode(fft_plan_forward_,
				CUFFT_COMPATIBILITY_NATIVE) != CUFFT_SUCCESS) {
    error->one(FLERR,"Failed to set compatibilit mode for forward FFT.");
  }
  if (cufftSetCompatibilityMode(fft_plan_inverse_,
				CUFFT_COMPATIBILITY_NATIVE) != CUFFT_SUCCESS) {
    error->one(FLERR,"Failed to set compatibilit mode for inverse FFT.");
  }

  linf_ = new double[nu_];
  if (cudaMalloc(&dev_linf_, nu_*sizeof(double)) != cudaSuccess) {
    error->one(FLERR,"cudaMalloc failed.");
  }
}


void transpose_row(int idim, int nx, int ny,
		   double **inbuf, cufftDoubleReal *outbuf)
{
  for (int i = 0; i < nx*ny; i++) {
    outbuf[i] = inbuf[i][idim];
  }
}


void transpose_row(int idim, int nx, int ny,
		   double_complex **inbuf, cufftDoubleComplex *outbuf)
{
  int k, nyy;

  nyy = ny/2+1;

  k = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nyy; j++) {
      outbuf[k].x = creal(inbuf[i*ny+j][idim]);
      outbuf[k].y = cimag(inbuf[i*ny+j][idim]);
      k++;
    }
  }
}


void buffer_full_to_buffer_half(int ndof, int nx, int ny,
                                double_complex **inbuf,
                                cufftDoubleComplex *outbuf)
{
  int ndof_sq = ndof*ndof;
  int nyy = ny/2+1;

  int k = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nyy; j++) {
      for (int n = 0; n < ndof_sq; n++) {
        outbuf[k].x = creal(inbuf[i*ny+j][n]);
        outbuf[k].y = cimag(inbuf[i*ny+j][n]);
        k++;
      }
    }
  }
}


void GFMDSolverStaticCUDA::set_kernel(StiffnessKernel *kernel,
				      bool normalize)
{
  if (screen && comm->me == 0)
    fprintf(screen, "Computing stiffness matrices...\n");

  /*
   * Get stiffness matrix
   */

  double_complex **operator_buffer;
  cufftDoubleComplex *complex_buffer;

  operator_buffer = create_complex_operator_buffer("operator_buffer");
  fill_phi_buffer(ndof, nx, xlo_loc, xhi_loc, ny, ylo_loc, yhi_loc,
		  kernel, operator_buffer, normalize, error);

#if 0
  for (int i = 0; i < nxy_loc; i++) {
    for (int x = 0; x < ndof; x++) {
      for (int y = 0; y < ndof; y++) {
	int j = x*ndof+y;
	//	if (x == 1 && y == 0 || x == 0 && y == 1 || x == 2 && y == 2 || x == 3 && y == 3 || x == 4 && y == 5 || x == 5 && y == 4)
	if (x == y)
	  operator_buffer[i][j] = -1.0/nxy_loc;
	else
	  operator_buffer[i][j] = 0.0;
      }
    }
  }
#endif
  
  /*
   * Upload operator buffer to device memory
   */
#if 1
  dev_q_operator_pitch_ = nbuffer_half_loc_;
  complex_buffer = (cufftDoubleComplex *)
    malloc(nbuffer_half_loc_*sizeof(cufftDoubleComplex));

  for (int idim = 0; idim < ndof_sq; idim++) {
    transpose_row(idim, nx, ny, operator_buffer, complex_buffer);

    if (cudaMemcpy(&dev_q_operator_[idim*dev_q_operator_pitch_],
		   complex_buffer,
		   nbuffer_half_loc_*sizeof(cufftDoubleComplex),
		   cudaMemcpyHostToDevice) != cudaSuccess) {
      char errstr[1024];

      sprintf(errstr, "cudaMemcpy to device failed: %s.",
	      cudaGetErrorString(cudaGetLastError()));
      error->one(FLERR,errstr);
    }
  }

  free(complex_buffer);
#else
  complex_buffer = (cufftDoubleComplex *)
    malloc(nbuffer_half_loc_*ndof_sq*sizeof(cufftDoubleComplex));

  buffer_full_to_buffer_half(ndof, nx, ny, operator_buffer, complex_buffer);

  CU_CHECK( cudaMemcpy(dev_q_operator_, complex_buffer,
                       nbuffer_half_loc_*ndof_sq*sizeof(cufftDoubleComplex),
                       cudaMemcpyHostToDevice) );

  free(complex_buffer);
#endif

  destroy_complex_operator_buffer(operator_buffer);

  /*
   * Get linear force contributions
   */

  if (nu_ > 1) {
    for (int i = 0; i < nu_; i++) {
      kernel->get_force_at_gamma_point(linf_);
    }

    if (cudaMemcpy(dev_linf_, linf_, nu_*sizeof(double),
		   cudaMemcpyHostToDevice) != cudaSuccess) {
      char errstr[1024];
    
      sprintf(errstr, "cudaMemcpy to device failed: %s.",
	      cudaGetErrorString(cudaGetLastError()));
      error->one(FLERR,errstr);
    }
  }

  if (screen && comm->me == 0)
    fprintf(screen, "...done\n");
}


double GFMDSolverStaticCUDA::post_force(void *dev_u_ptr, void *dev_f_ptr,
					char *dump_prefix)
{
  double epot = 0.0;

  cufftDoubleReal *dev_u = static_cast<cufftDoubleReal*>(dev_u_ptr);
  cufftDoubleReal *dev_f = static_cast<cufftDoubleReal*>(dev_f_ptr);

  size_t dev_u_pitch = nxy_loc;
  size_t dev_f_pitch = nxy_loc;

  if (dev_u == dev_f) {
    error->all(FLERR,"dev_u == dev_f. In-place operation not supported.");
  }

  if (!cu_on_device(dev_u)) {
    error->one(FLERR,"dev_u needs to be on device.");
  }
  if (!cu_on_device(dev_f)) {
    error->one(FLERR,"dev_f needs to be on device.");
  }

  /*
   * Perform forward FFT for each dimension
   */
#ifdef AVOID_BATCH_TRANSFORM
  for (int idim = 0; idim < ndof; idim++) {
    if (cufftExecD2Z(fft_plan_forward_,
		     dev_u          + idim*dev_u_pitch,
		     dev_q_buffer2_ + idim*dev_q_buffer_pitch_
		     ) != CUFFT_SUCCESS) {
      error->all(FLERR,"Forward FFT failed.");
    }
  }
#else
  if (cufftExecD2Z(fft_plan_forward_, dev_u, dev_q_buffer2_) != CUFFT_SUCCESS) {
    error->all(FLERR,"Forward FFT failed.");
  }
#endif
  
#if 0
  // just copy input to output buffer --- for debugging purposes
  cudaMemcpy(dev_q_buffer1_, dev_q_buffer2_,
	     dev_q_buffer_pitch_*ndof*sizeof(cufftDoubleComplex),
	     cudaMemcpyDeviceToDevice);
#endif

  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
   // dev_q_buffer1 is f (output) and dev_q_buffer2 is u (input)
  if (gfmd_solver_cuda_apply_operator(nbuffer_half_loc_, ndof,
				      dev_q_operator_, dev_q_operator_pitch_,
				      dev_q_buffer2_, dev_q_buffer1_,
				      dev_q_buffer_pitch_, gammai_, dev_linf_)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "gfmd_solver_cuda_kernel failed: %s.",
	    cudaGetErrorString(cudaGetLastError()));
    error->all(FLERR,errstr);
  }

  /*
   * Add gamma point energy. 1/2 of this energy is already included below.
   */
  if (gammai_ >= 0 && nu_ > 1) {
    for (int i = 0; i < nu_; i++) {
      double f0;
      cudaMemcpy(&f0, &dev_q_buffer2_[(3*i+2)*dev_q_buffer_pitch_+gammai_],
		 sizeof(double), cudaMemcpyDeviceToHost);
      epot -= 0.5*linf_[i]*f0; //creal(q_buffer_[gammai_][3*i+2]);
    }
  }

  /*
   * Perform inverse FFT for each dimension
   */
#ifdef AVOID_BATCH_TRANSFORM
  for (int idim = 0; idim < ndof; idim++) {
    if (cufftExecZ2D(fft_plan_inverse_,
		     dev_q_buffer1_ + idim*dev_q_buffer_pitch_,
		     dev_f          + idim*dev_f_pitch
		     ) != CUFFT_SUCCESS) {
      error->all(FLERR,"Inverse FFT failed.");
    }
  }
#else
  if (cufftExecZ2D(fft_plan_inverse_, dev_q_buffer1_, dev_f) != CUFFT_SUCCESS) {
    error->all(FLERR,"Inverse FFT failed.");
  }
#endif

  /*
   * Compute potential energy
   */
  if (gfmd_solver_cuda_epot(nxy_loc, ndof,
			    dev_u, dev_u_pitch,
			    dev_f, dev_f_pitch,
			    dev_epot_buf_, epot)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "gfmd_grid_cuda_epot failed: %s.",
	    cudaGetErrorString(cudaGetLastError()));
    error->all(FLERR,errstr);
  }

  return epot;
}
