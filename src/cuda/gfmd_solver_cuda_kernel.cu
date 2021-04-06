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

#define THREADS_PER_BLOCK 256

#define CU_CHECK(x)  { cudaError err = x; if (err != cudaSuccess) return err; }

#define CMUL_REAL(a, b)  ((a).x*(b).x - (a).y*(b).y)
#define CMUL_IMAG(a, b)  ((a).x*(b).y + (a).y*(b).x)


template<typename T>
__device__ void warp_reduce(volatile T *sdata, int tid)
{
  sdata[tid] += sdata[tid + 32];
  sdata[tid] += sdata[tid + 16];
  sdata[tid] += sdata[tid +  8];
  sdata[tid] += sdata[tid +  4];
  sdata[tid] += sdata[tid +  2];
  sdata[tid] += sdata[tid +  1];
}


/*
 * GFMD CUDA kernel: Multiply operator to device buffer for each q-vector.
 * Fixed number of degrees of freedom
 */
template<unsigned int blockDim_x, unsigned int ndof>
__global__ void mul_op_buf(int nbuf, cufftDoubleComplex *op,
                           int op_pitch, cufftDoubleComplex *inbuf,
                           cufftDoubleComplex *outbuf, int buf_pitch)
{
  int q_i   = (blockIdx.x / ndof)*blockDim_x + threadIdx.x;
  int dof_i =  blockIdx.x % ndof;

  //int i = blockIdx.x*blockDim_x + threadIdx.x;
  //int q_i = i / ndof;
  //int dof_i = i % ndof;

  if (q_i < nbuf && dof_i < ndof) {
    int op_i = q_i + dof_i*ndof*op_pitch;
    //int op_i = (q_i*ndof + dof_i)*ndof;

    cufftDoubleReal v_real = 0.0;
    cufftDoubleReal v_imag = 0.0;
    cufftDoubleComplex v_op;
    cufftDoubleComplex v_buf;

    //if (q_i == 4) printf("id4: phi=%f  u=%f  \n",op   [op_i + 0 ],inbuf[ q_i + 0]);//TAS
    // force = stiffness x displacement
    for (int dof_j = 0; dof_j < ndof; dof_j++) {
      v_op  = op   [op_i + dof_j*op_pitch ];
      //cufftDoubleComplex v_op  = op   [op_i + dof_j          ];
      v_buf = inbuf[ q_i + dof_j*buf_pitch];
    


      v_real -= CMUL_REAL(v_op, v_buf);
      v_imag -= CMUL_IMAG(v_op, v_buf);
    }
    //if (q_i == 4) printf("id4: phi=%f  u=%f  f=%f\n",op   [op_i + 0 ],inbuf[ q_i + 0]);//TAS

    outbuf[q_i + dof_i*buf_pitch].x = v_real;
    outbuf[q_i + dof_i*buf_pitch].y = v_imag;
  }
}


/*
 * GFMD CUDA kernel: Multiply operator to device buffer for each q-vector
 */
template<unsigned int blockDim_x>
__global__ void mul_op_buf_gen(int nbuf, int ndof, cufftDoubleComplex *op,
                               int op_pitch, cufftDoubleComplex *inbuf,
                               cufftDoubleComplex *outbuf, int buf_pitch)
{
  int q_i   = (blockIdx.x / ndof)*blockDim_x + threadIdx.x;
  int dof_i =  blockIdx.x % ndof;

  if (q_i < nbuf && dof_i < ndof) {
    int op_i = q_i + dof_i*ndof*op_pitch;
    //int op_i = q_i*op_pitch + dof_i*ndof;

    cufftDoubleReal v_real = 0.0;
    cufftDoubleReal v_imag = 0.0;

    for (int dof_j = 0; dof_j < ndof; dof_j++) {
      cufftDoubleComplex v_op  = op   [op_i + dof_j*op_pitch ];
      //cufftDoubleComplex v_op  = op   [op_i + dof_j          ];
      cufftDoubleComplex v_buf = inbuf[ q_i + dof_j*buf_pitch];

      v_real -= CMUL_REAL(v_op, v_buf);
      v_imag -= CMUL_IMAG(v_op, v_buf);
    }

    outbuf[q_i + dof_i*buf_pitch].x = v_real;
    outbuf[q_i + dof_i*buf_pitch].y = v_imag;
  }
}


/*
 * GFMD CUDA kernel: Add linear force contribution to q=0 component
 */
__global__ void add_linf(int nu, double *linf, int gammai, 
			 cufftDoubleComplex *f, int buf_pitch)
{
  int i = threadIdx.x;

  f[(3*i+2)*buf_pitch+gammai].x += linf[i];
}


__global__ void reduce_mul(int nbuf, cufftDoubleReal *buf1,
			   cufftDoubleReal *buf2, cufftDoubleReal *bufout)
{
  extern __shared__ cufftDoubleReal sdata[];

  int tid = threadIdx.x;
  int i = 2*blockIdx.x*blockDim.x + tid;

  if (i+blockDim.x < nbuf) {
    sdata[tid] = buf1[i]*buf2[i] + buf1[i+blockDim.x]*buf2[i+blockDim.x];
  }
  else if (i < nbuf) {
    sdata[tid] = buf1[i]*buf2[i];
  }
  else {
    sdata[tid] = 0.0;
  }
  __syncthreads();

  for (int s = blockDim.x/2; s > 32; s >>= 1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  if (tid < 32)
    warp_reduce(sdata, tid);

  if (tid == 0)
    bufout[blockIdx.x] = sdata[0];
}

