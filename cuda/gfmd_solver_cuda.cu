
#include "gfmd_solver_cuda_cu.h"
#include "gfmd_solver_cuda_kernel.cu"

#define THREADS_PER_BLOCK 256

#define CU_CHECK(x)  { cudaError err = x; if (err != cudaSuccess) return err; }

#define CMUL_REAL(a, b)  ((a).x*(b).x - (a).y*(b).y)
#define CMUL_IMAG(a, b)  ((a).x*(b).y + (a).y*(b).x)


/*
 * GFMD CUDA kernel interface
 */
cudaError gfmd_solver_cuda_apply_operator(int nbuf, int ndof,
					  cufftDoubleComplex *op,
					  int op_pitch,
					  cufftDoubleComplex *u,
					  cufftDoubleComplex *f,
					  int buf_pitch, int gammai, 
					  double *linf)
{
  int num_blocks = nbuf*ndof / THREADS_PER_BLOCK;

  while (num_blocks*THREADS_PER_BLOCK < nbuf*ndof)
    num_blocks++;
  if (num_blocks % ndof != 0)
    num_blocks = ((num_blocks/ndof)+1)*ndof;

  /*
   * Multiply stiffness to displacements to get the forces
   */
  switch (ndof) {
  case 3:
    mul_op_buf<THREADS_PER_BLOCK, 3><<<num_blocks, THREADS_PER_BLOCK>>>
      (nbuf, op, op_pitch, u, f, buf_pitch);
    break;
  case 6:
    mul_op_buf<THREADS_PER_BLOCK, 6><<<num_blocks, THREADS_PER_BLOCK>>>
      (nbuf, op, op_pitch, u, f, buf_pitch);
    break;
  case 9:
    mul_op_buf<THREADS_PER_BLOCK, 9><<<num_blocks, THREADS_PER_BLOCK>>>
      (nbuf, op, op_pitch, u, f, buf_pitch);
    break;
  default:
    mul_op_buf_gen<THREADS_PER_BLOCK><<<num_blocks, THREADS_PER_BLOCK>>>
      (nbuf, ndof, op, op_pitch, u, f, buf_pitch);
    break;
  }
  CU_CHECK( cudaPeekAtLastError() );

  /*
   * Add linear force contribution to q=0 components
   */
  int nu = ndof/3;
  if (gammai >= 0 && nu > 1) {//printf("in the kernel adding linf[0]:\n");//TAS
    add_linf<<<1, nu>>>(nu, linf, gammai, f, buf_pitch);
    CU_CHECK( cudaPeekAtLastError() );
  }

  return cudaSuccess;
}


/*
 * Compute potential energy
 * Note: u will be overriden
 */
cudaError gfmd_solver_cuda_epot(int nbuf, int ndof,
				cufftDoubleReal *dev_u, int dev_u_pitch,
				cufftDoubleReal *dev_f, int dev_f_pitch,
				cufftDoubleReal *&dev_epot_buf,
				double &epot_out)
{
  int idim, j;
  int num_blocks = nbuf / (2*THREADS_PER_BLOCK);

  cufftDoubleReal epot_buf[nbuf];
  double epot = 0.0;

#if 0
  cufftDoubleReal _u[nbuf*ndof], _f[nbuf*ndof];
  cudaMemcpy(_u, dev_u, nbuf*ndof*sizeof(cufftDoubleReal),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(_f, dev_f, nbuf*ndof*sizeof(cufftDoubleReal),
	     cudaMemcpyDeviceToHost);
  epot = 0.0;
  for (j = 0; j < nbuf*ndof; j++)
    epot += _u[j]*_f[j];
  epot *= -0.5;
  printf("epot = %f\n", epot);

  printf("2: u - %e %e %e %e %e %e\n", _u[0*nbuf+10], _u[1*nbuf+20],
	 _u[2*nbuf+30], _u[3*nbuf+40], _u[4*nbuf+50], _u[5*nbuf+60]);
  printf("2: f - %e %e %e %e %e %e\n", _f[0*nbuf+10], _f[1*nbuf+20],
	 _f[2*nbuf+30], _f[3*nbuf+40], _f[4*nbuf+50], _f[5*nbuf+60]);
#endif

  if (num_blocks*2*THREADS_PER_BLOCK < nbuf)
    num_blocks++;

  if (!dev_epot_buf) {
    CU_CHECK( cudaMalloc(&dev_epot_buf, num_blocks*sizeof(cufftDoubleReal)) );//printf("#CUDA copying memory to host gfmdsolver_cuda_epot\n");//TAS
  }

  /*
   * Call reduction kernel and sum potential energy
   */
  for (idim = 0; idim < ndof; idim++) {
    reduce_mul
      <<<num_blocks, THREADS_PER_BLOCK,
      THREADS_PER_BLOCK*sizeof(cufftDoubleReal)>>>
      (nbuf, dev_u + idim*dev_u_pitch, dev_f + idim*dev_f_pitch, dev_epot_buf);
    CU_CHECK( cudaPeekAtLastError() );

    CU_CHECK( cudaMemcpy(epot_buf, dev_epot_buf,
			 num_blocks*sizeof(cufftDoubleReal),
			 cudaMemcpyDeviceToHost) );//printf("#CUDA copying memory to host gfmdsolver_cuda_epot\n");//TAS

    for (j = 0; j < num_blocks; j++) {
      epot += epot_buf[j];
    }
  }

  epot_out -= 0.5*epot;
  return cudaSuccess;
}
