#ifdef GFMD_CUFFT
#include <pthread.h>
#endif

#include <string.h>
#include <stdlib.h>

#include "comm.h"
#include "domain.h"
#include "lammps.h"
#include "memory.h"
#include "mpi.h"
#include "timer.h"

#include "linearalgebra.h"

#include "gfmd_memory.h"
#include "gfmd_solver.h"

#ifdef GFMD_CUFFT
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>

#include "gfmd_solver_cuda.h"
#endif

using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))


/* ----------------------------------------------------------------------
 * Base class for 2d GFMD grids
 * --------------------------------------------------------------------*/

GFMDSolver::GFMDSolver(LAMMPS *lmp, size_t in_nx, size_t in_ny, size_t in_ndof)
  : Pointers(lmp)
{
  strcpy(name, "none");

  nx = in_nx;
  ny = in_ny;
  ndof = in_ndof;

  ndof_sq = ndof*ndof;

  /*
   * Get basic MPI information
   */
  me = comm->me;
  nprocs = comm->nprocs;

  /*
   * Check which part of the grid should be stored on this processor
   */
  xlo_loc = (size_t) (nx*domain->sublo[0]/domain->xprd);
  xhi_loc = (size_t) (nx*domain->subhi[0]/domain->xprd-1);

  ylo_loc = (size_t) (ny*domain->sublo[1]/domain->yprd);
  yhi_loc = (size_t) (ny*domain->subhi[1]/domain->yprd-1);

  nx_loc  = xhi_loc-xlo_loc+1;
  ny_loc  = yhi_loc-ylo_loc+1;

  nxy_loc = nx_loc*ny_loc;

  //  printf("me = %i, xlo_loc = %zu, xhi_loc = %zu, ylo_loc = %zu, "
  //	 "xhi_loc = %zu\n", me, xlo_loc, xhi_loc, ylo_loc, yhi_loc);
}


GFMDSolver::~GFMDSolver()
{
}


double **GFMDSolver::create_double_buffer(char *name)
{
  return memory->create_2d_double_array(MAX(1, nxy_loc), ndof, name);
}


void GFMDSolver::destroy_double_buffer(double **buffer)
{
  memory->destroy_2d_double_array(buffer);
}


double_complex **GFMDSolver::create_complex_buffer(char *name)
{
  return create_2d_complex_array(memory, MAX(1, nxy_loc), ndof, name);
}


void GFMDSolver::destroy_complex_buffer(double_complex **buffer)
{
  destroy_2d_complex_array(memory, buffer);
}


double_complex **GFMDSolver::create_complex_operator_buffer(char *name)
{
  return create_2d_complex_array(memory, MAX(1, nxy_loc), ndof_sq, name);
}


void GFMDSolver::destroy_complex_operator_buffer(double_complex **buffer)
{
  destroy_2d_complex_array(memory, buffer);
}


double GFMDSolver::memory_usage()
{
  return 0.0;
}


/* ----------------------------------------------------------------------
 * Default implementation using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

GFMDSolverDefault::GFMDSolverDefault(LAMMPS *lmp,
				     size_t in_nx, size_t in_ny,
				     size_t in_ndof)
  : GFMDSolver(lmp, in_nx, in_ny, in_ndof)
{
  int nfft;

  strcpy(name, "default");

  fft  = new FFT3d(lmp, world,
	       //  fast, med,              slow
		   1,    ny,               nx,
		   0, 0, ylo_loc, yhi_loc, xlo_loc, xhi_loc,
		   0, 0, ylo_loc, yhi_loc, xlo_loc, xhi_loc,
		   0, 0, &nfft);

  /*
  if (nfft != nxy_loc) {
    char errstr[1024];
    sprintf(errstr, "Local buffer size (= %zu) does not match buffer size "
	    "return by the FFT wrapper (= %i).", nxy_loc, nfft);
    error->one(errstr);
  }
  */

  fft_data = (double *) memory->smalloc(nxy_loc*2*sizeof(double),
					"GFMDSolver::fft_data");

  q_buffer = create_complex_buffer("GFMDSolver::q_buffer");
  q_operator = NULL;
}


GFMDSolverDefault::~GFMDSolverDefault()
{
  destroy_complex_buffer(q_buffer);

  memory->sfree(fft_data);

  delete fft;
}


void GFMDSolverDefault::set_operator(double_complex **operator_buffer)
{
  q_operator = operator_buffer;
}


double GFMDSolverDefault::join(double **input_buffer,
					 double **output_buffer)
{
  double epot, tmp;

  double_complex F_q[ndof];

  /*
   * Fill FFT buffer and perform FFT on local data
   */
  for (int idim = 0; idim < ndof; idim++) {
    int m = 0;
    for (int idx = 0; idx < nxy_loc; idx++) {
      /*
       * displacement U is all real.  space it out into everyother element
       */
      fft_data[m++] = input_buffer[idim][idx];
      fft_data[m++] = 0.;
    }
    /*
     * perform the FFT!
     */
    fft->compute(fft_data, fft_data, -1);
    m = 0;
    for (int idq = 0; idq < nxy_loc; idq++){
      // pack double array data into complex array
      q_buffer[idq][idim] = fft_data[m] + fft_data[m+1]*I;
      m += 2;
    }
  }

  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
  epot = 0.0;
  for (int idq = 0; idq < nxy_loc; idq++) {
    MatMulVec(ndof, q_operator[idq], q_buffer[idq], F_q);
    for (int idim = 0; idim < ndof; idim++) {
      epot += creal(conj(F_q[idim]) * q_buffer[idq][idim]);
      q_buffer[idq][idim] = -F_q[idim];
    }
  }
  epot *= 0.5;

  /*
   * Gather total energy on all processors
   */
  MPI_Allreduce(&epot, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  epot = tmp;

  /*
   * Now transform F(q) to F(r)
   */
  for (int idim = 0; idim < ndof; idim++){
    int m = 0;
    for (int idq = 0; idq < nxy_loc; idq++){
      fft_data[m++] = creal(q_buffer[idq][idim]);
      fft_data[m++] = cimag(q_buffer[idq][idim]);
    }
    fft->compute(fft_data, fft_data, 1);
    m = 0;
    //    for (int idx = 0; idx < nxy_loc; idx++){
    for (int ix = xlo_loc; ix <= xhi_loc; ix++) {
      for (int iy = ylo_loc; iy <= yhi_loc; iy++) {
	int idx = ix*ny+iy;
	output_buffer[idim][idx] = fft_data[m];
	m += 2;
      }
    }
      //    }
  }

  return epot;
}


double GFMDSolverDefault::memory_usage()
{
  double bytes = 0.0;

  // fft_data
  bytes += nxy_loc*2*sizeof(double);

  return bytes;
}


/* ----------------------------------------------------------------------
 * Use GFMD_CUFFT
 * --------------------------------------------------------------------*/

#ifdef GFMD_CUFFT
// all cuda stuff should be handled in a separate thread, such that CPU
// computations can continue will the GPU is working
// the mutex is used for synchronization
CUDAThread::CUDAThread(LAMMPS *lmp, size_t in_nx, size_t in_ny, size_t in_ndof,
		       double_complex **operator_buffer)
  : GFMDSolver(lmp, in_nx, in_ny, in_ndof)
{
  char *cuda_device_str, *endptr;

  q_operator = operator_buffer;

  cuda_device_id = 0;
  cuda_device_str = getenv("CUDA_DEVICE");
  if (cuda_device_str) {
    cuda_device_id = strtol(cuda_device_str, &endptr, 10);
    if (endptr == cuda_device_str) {
      char errstr[1024];
      sprintf(errstr, "[CUDAThread::CUDAThread] Found CUDA_DEVICE=%s, "
	      "but could not convert to integer device id.", cuda_device_str);
      error->one(errstr);
    }
    if (screen)
      fprintf(screen, "Using CUDA device # %i\n", cuda_device_id);
  }

  if (pthread_mutex_init(&mutex, NULL)) {
    error->one("[CUDAThread::CUDAThread] Failed to create mutex.");
  }
  if (pthread_cond_init(&cond, NULL)) {
    error->one("[CUDAThread::CUDAThread] Failed to create condition.");
  }

  thread_is_running = true;
}


CUDAThread::~CUDAThread()
{
  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond);

  cufftDestroy(fft_plan_forward);
  cufftDestroy(fft_plan_inverse);

  cudaFree(dev_q_operator);
  cudaFree(dev_u);
  cudaFree(dev_q_buffer1);
  cudaFree(dev_q_buffer2);

  if (dev_epot_buf)
    cudaFree(dev_epot_buf);

  free(complex_buffer);
}


void CUDAThread::init()
{
  size_t pitch;

  if (cudaSetDevice(cuda_device_id) != cudaSuccess) {
    error->all("[CUDAThread::CUDAThread] cudaSetDevice failed.");
  }

  dev_epot_buf = NULL;

  nbuffer_half_loc = nx*(ny/2+1);

  if (nprocs != 1) {
    error->all("[CUDAThread::CUDAThread] CUDA FFT for GFMD "
	       "currently only works on a single processor");
  }

  /* Don't know if this is always true */
  if (sizeof(cufftDoubleReal) != sizeof(double)) {
    error->all("[CUDAThread::CUDAThread] cufftDoubleReal and double "
	       "size mismatch.");
  }
  if (sizeof(cufftDoubleComplex) != sizeof(double_complex)) {
    error->all("[CUDAThread::CUDAThread] cufftDoubleComplex and "
	       "double_complex size mismatch.");
  }

  if (cudaMallocPitch(&dev_q_operator, &dev_q_operator_pitch,
		      nbuffer_half_loc*sizeof(cufftDoubleComplex), ndof_sq)
      != cudaSuccess) {
    error->all("[CUDAThread::CUDAThread] cudaMalloc failed.");
  }

  if (cudaMallocPitch(&dev_q_buffer1, &dev_q_buffer_pitch,
		      nbuffer_half_loc*sizeof(cufftDoubleComplex), ndof)
      != cudaSuccess) {
    error->all("[CUDAThread::CUDAThread] cudaMalloc failed.");
  }

  if (cudaMallocPitch(&dev_q_buffer2, &pitch,
		      nbuffer_half_loc*sizeof(cufftDoubleComplex), ndof)
      != cudaSuccess) {
    error->all("[CUDAThread::CUDAThread] cudaMalloc failed.");
  }

  if (pitch != dev_q_buffer_pitch) {
    error->all("[CUDAThread::CUDAThread] pitch != dev_q_buffer_pitch");
  }

  pitch = dev_q_operator_pitch/sizeof(cufftDoubleComplex);
  if (pitch*sizeof(cufftDoubleComplex) != dev_q_operator_pitch) {
    error->all("[CUDAThread::CUDAThread] pitch*sizeof("
	       "cufftDoubleComplex) != dev_q_operator_pitch");
  }
  dev_q_operator_pitch = pitch;

  pitch = dev_q_buffer_pitch/sizeof(cufftDoubleComplex);
  if (pitch*sizeof(cufftDoubleComplex) != dev_q_buffer_pitch) {
    error->all("[CUDAThread::CUDAThread] pitch*sizeof("
	       "cufftDoubleComplex) != dev_q_buffer_pitch");
  }
  dev_q_buffer_pitch = pitch;

  //  dev_real_buffer1 = (cufftDoubleReal *) dev_q_buffer1;
  dev_f = (cufftDoubleReal *) dev_q_buffer2;
  dev_f_pitch = 2*dev_q_buffer_pitch;

  if (cudaMallocPitch(&dev_u, &dev_u_pitch,
		      nxy_loc*sizeof(cufftDoubleReal), ndof)
      != cudaSuccess) {
    error->all("[CUDAThread::CUDAThread] cudaMalloc failed.");
  }

  pitch = dev_u_pitch/sizeof(cufftDoubleReal);
  if (pitch*sizeof(cufftDoubleReal) != dev_u_pitch) {
    error->all("[CUDAThread::CUDAThread] pitch*sizeof("
	       "cufftDoubleComplex) != dev_u_pitch");
  }
  dev_u_pitch = pitch;

  //  printf("pitches: %i %i\n", dev_q_operator_pitch, dev_q_buffer_pitch);

  complex_buffer = (cufftDoubleComplex *)
    malloc(dev_q_buffer_pitch*ndof*sizeof(cufftDoubleComplex));
  real_buffer = (cufftDoubleReal *) complex_buffer;

  if (cufftPlan2d(&fft_plan_forward, nx, ny, CUFFT_D2Z) != CUFFT_SUCCESS) {
    error->all("[CUDAThread::CUDAThread] Failed to plan forward FFT.");
  }
  if (cufftPlan2d(&fft_plan_inverse, nx, ny, CUFFT_Z2D) != CUFFT_SUCCESS) {
    error->all("[CUDAThread::CUDAThread] Failed to plan inverse FFT.");
  }

  set_operator(q_operator);
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


void CUDAThread::set_operator(double_complex **operator_buffer)
{
  /*
   * Upload operator buffer to device memory
   */
  for (int idim = 0; idim < ndof_sq; idim++) {
    transpose_row(idim, nx, ny, operator_buffer, complex_buffer);

    if (cudaMemcpy(&dev_q_operator[idim*dev_q_operator_pitch], complex_buffer,
		   nbuffer_half_loc*sizeof(cufftDoubleComplex),
		   cudaMemcpyHostToDevice) != cudaSuccess) {
      char errstr[1024];

      sprintf(errstr, "[GFMDSolverCUDA::set_operator] cudaMemcpy to device "
	      "failed: %s.", cudaGetErrorString(cudaGetLastError()));
    }
  }
}


double CUDAThread::apply_operator()
{
  int i, j;
  double epot;

  double_complex F_q[ndof];

  if (input_buffer == output_buffer) {
    error->all("[GFMDSolverCUDA::join] input_buffer == "
	       "output_buffer");
  }

  /*
   * Copy data to device
   */
  for (int idim = 0; idim < ndof; idim++) {
    if (cudaMemcpy(dev_u + idim*dev_u_pitch, input_buffer[idim],
		   nxy_loc*sizeof(cufftDoubleReal),
		   cudaMemcpyHostToDevice) != cudaSuccess) {
      char errstr[1024];

      sprintf(errstr, "[GFMDSolverCUDA::join] cudaMemcpy to device "
	      "failed: %s.", cudaGetErrorString(cudaGetLastError()));
      error->all(errstr);
    }
  }

  /*
   * Perform forward FFT for each dimension
   */
  for (int idim = 0; idim < ndof; idim++) {
    if (cufftExecD2Z(fft_plan_forward,
		     dev_u         + idim*dev_u_pitch,
		     dev_q_buffer2 + idim*dev_q_buffer_pitch
		     ) != CUFFT_SUCCESS) {
      error->all("[GFMDSolverCUDA::join] Forward FFT failed.");
    }
  }

  /*
   * Perform matrix operation: F(q) = -Phi(q) x U(q)
   * and compute potential energy
   */
  if (gfmd_solver_cuda_apply_operator(nbuffer_half_loc, ndof,
				      dev_q_operator, dev_q_operator_pitch,
				      dev_q_buffer2, dev_q_buffer1,
				      dev_q_buffer_pitch)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "[GFMDSolverCUDA::join] gfmd_solver_cuda_kernel "
	    "failed: %s.", cudaGetErrorString(cudaGetLastError()));
    error->all(errstr);
  }
		 
  /*
   * Perform inverse FFT for each dimension
   */
  for (int idim = 0; idim < ndof; idim++) {
    if (cufftExecZ2D(fft_plan_inverse,
		     dev_q_buffer1 + idim*dev_q_buffer_pitch,
		     dev_f         + idim*dev_f_pitch
		     ) != CUFFT_SUCCESS) {
      error->all("[GFMDSolverCUDA::join] Inverse FFT failed.");
    }
  }

  /*
   * Compute potential energy
   */
  if (gfmd_solver_cuda_epot(nxy_loc, ndof,
			    dev_u, dev_u_pitch,
			    dev_f, dev_f_pitch,
			    &dev_epot_buf, &epot)
      != cudaSuccess) {
    char errstr[1024];

    sprintf(errstr, "[GFMDSolverCUDA::join] gfmd_grid_cuda_epot "
	    "failed: %s.", cudaGetErrorString(cudaGetLastError()));
    error->all(errstr);
  }

  /*
   * Copy data from device
   */
  for (int idim = 0; idim < ndof; idim++) {
    if (cudaMemcpy(real_buffer,
		   dev_f + idim*dev_f_pitch,
		   nxy_loc*sizeof(cufftDoubleReal),
		   cudaMemcpyDeviceToHost) != cudaSuccess) {
      char errstr[1024];

      sprintf(errstr, "[GFMDSolverCUDA::join] cudaMemcpy from "
	      "device failed: %s.",
	      cudaGetErrorString(cudaGetLastError()));
      error->all(errstr);
    }

    for (int ix = 0; ix < nx_loc; ix++) {
      memcpy(&output_buffer[idim][(ix+xlo_loc)*ny + ylo_loc],
	     &real_buffer[ix*ny_loc],
	     nx_loc*sizeof(double));
    }
  }

  return epot;
}


void CUDAThread::run()
{
  init();

  while (thread_is_running) {
    lock();
    wait();
    if (thread_is_running) {
      epot = apply_operator();
    }
    unlock();
  }
}


void *GFMDSolverCUDA::cuda_thread_dispatcher(void *self_ptr)
{
  CUDAThread *self = static_cast<CUDAThread *>(self_ptr);

  self->run();

  delete self;
}
#endif


GFMDSolverCUDA::GFMDSolverCUDA(LAMMPS *lmp,
			       size_t in_nx, size_t in_ny, size_t in_ndof)
  : GFMDSolver(lmp, in_nx, in_ny, in_ndof)
{
#ifdef GFMD_CUFFT
  strcpy(name, "cuda");

  thread = NULL;
#else
  error->all("CUDA FFT for GFMD not available. Compile with -DGFMD_CUFFT.");
#endif
}


GFMDSolverCUDA::~GFMDSolverCUDA()
{
#ifdef GFMD_CUFFT
  thread->stop();
#endif
}


void GFMDSolverCUDA::detach(double **input_buffer, double **output_buffer)
{
#ifdef GFMD_CUFFT
  thread->input_buffer = input_buffer;
  thread->output_buffer = output_buffer;

  // tell thread to start solving
  thread->signal();
  thread->unlock();
#endif
}


double GFMDSolverCUDA::join(double **input_buffer, double **output_buffer)
{
#ifdef GFMD_CUFFT
  double tmp;

  // wait for thread to finish
  thread->lock();

  // Gather total energy on all processors
  MPI_Allreduce(&thread->epot, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);

  return tmp;
#else
  return 0.0;
#endif
}


void GFMDSolverCUDA::set_operator(double_complex **operator_buffer)
{
#ifdef GFMD_CUFFT
  if (thread) {
    error->all("GFMDSolverCUDA::set_operator may only be called once.");
  }

  thread = new CUDAThread(lmp, nx, ny, ndof, operator_buffer);
  thread->lock();
  if (pthread_create(&pthread, NULL,
		     GFMDSolverCUDA::cuda_thread_dispatcher, thread)) {
    error->one("GFMDSolverCUDA::GFMDSolverCUDA: Failed to create thread.");
  }

#else
  error->all("CUDA FFT for GFMD not available. Compile with -DGFMD_CUFFT.");
#endif
}



namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Instantiate a GFMD grid according to keyword arguments
 * --------------------------------------------------------------------*/

GFMDSolver *gfmd_solver_factory(char *keyword, LAMMPS *lmp,
				size_t nx, size_t ny, size_t ndof)
{
  GFMDSolver *grid = NULL;

  if (keyword) {
    if (!strcmp(keyword, "default")) {
      grid = new GFMDSolverDefault(lmp, nx, ny, ndof);
    }
    else if (!strcmp(keyword, "cuda")) {
      grid = new GFMDSolverCUDA(lmp, nx, ny, ndof);
    }
  }
  else {
    grid = new GFMDSolverDefault(lmp, nx, ny, ndof);
  }

  if (grid && keyword) {
    if (strcmp(grid->get_name(), keyword)) {
      char errstr[120];

      sprintf(errstr, "gfmd_solver_factory: Internal error: keyword '%s' "
	      "and fft grid class '%s' name mismatch.", keyword,
	      grid->get_name());
      lmp->error->all(errstr);
    }
  }

  return grid;
}

}
