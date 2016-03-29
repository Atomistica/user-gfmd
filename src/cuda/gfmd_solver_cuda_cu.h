#ifndef GFMD_SOLVER_CUDA_CU_H
#define GFMD_SOLVER_CUDA_CU_H

#include <cufft.h>

extern "C" {

cudaError gfmd_solver_cuda_apply_operator(int, int, cufftDoubleComplex *,
					  int, cufftDoubleComplex *,
					  cufftDoubleComplex *, int, int,
					  double *);

cudaError gfmd_solver_cuda_epot(int, int, cufftDoubleReal *, int,
				cufftDoubleReal *, int, cufftDoubleReal *&,
				double &);

};

#endif
