#ifndef FIX_FIRE_CUDA_CU_H
#define FIX_FIRE_CUDA_CU_H

//#include <cuda_runtime.h>
//
//cudaError fix_fire_cuda_dot_products(int, int, int, int *, double *, double *,
//				     double &, double &, double &);
//
//cudaError fix_fire_cuda_mix(double, double, int, int, int, int *, double *,
//			    double *);
#include "cuda_shared.h"

//extern "C" template<typename T>
//cudaError reduce(int n, T *arr, T &res);

//#include "cuda_shared.h"     // nve has it
//#include "cuda_common.h"
//#include "crm_cuda_utils.cu" // nve has it

//#include "fix_fire_cuda_cu.h"
//#include "fix_fire_cuda_kernel.cu"


/*!
 * Launcher for array reduction
 */
template<typename T>
cudaError reduce(int n, T *arr, T &res);

/*!
 * Compute the dot products v.f, v.v and f.f
 */
cudaError fix_fire_cuda_dot_products(int groupbit, int nlocal, int nmax,
                                     int *dev_mask, double *dev_v,
                                     double *dev_f, double &vf,
                                     double &vg_dot_vg, double &Fg_dot_Fg);

/*!
 * Mix each atom's force vector into its velocity vector
 */
cudaError fix_fire_cuda_mix(double a, double b, int groupbit, int nlocal,
                            int nmax, int *dev_mask, double *dev_v,
                            double *dev_f);

#endif
