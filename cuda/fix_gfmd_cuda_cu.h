#ifndef FIX_GFMD_CUDA_CU_H
#define FIX_GFMD_CUDA_CU_H

#include <cuda_runtime.h>
//#include <mpi.h>

#include "../../src/lmptype.h"

using namespace LAMMPS_NS;

extern "C" {

cudaError fix_gfmd_cuda_list_to_grid(double xprd, double yprd, 
				     double xprd_half, double yprd_half,
				     int groupbit, int nlocal, int nall,
				     int nmax, int *mask, X_FLOAT *x, int *gid,
				     X_FLOAT *xeq, int dxshift, int dyshift,
				     int nx, int ny, int xlo_loc, int ylo_loc,
				     int nx_loc, int ny_loc, int nxy_loc,
				     double *u_xy, int &natoms_cur,
				     int &ngfmd_loc);

cudaError fix_gfmd_cuda_grid_to_list(int groupbit, int nlocal, int nall,
				     int nmax, int *dev_mask, F_FLOAT *dev_data,
				     int *dev_gid, int nx, int ny, int nxy,
				     int xlo_loc, int ylo_loc, int xhi_loc,
				     int yhi_loc, double *dev_grid,
				     bool *dev_force_is_known, int &natoms_cur,
				     double sum_loc[3]);

cudaError fix_gfmd_cuda_xcm(int groupbit, double xprd,
			    double yprd, double zprd, int nlocal, int nmax,
			    int *mask, int *type, double *mass, double *rmass,
			    double *x, tagint *image,
			    double cmone[3]);

cudaError CudaData_printmyaddress(void* device_ptr);// TAS
//void CudaData_printmyaddre(float* device_ptr);// TAS
//void CudaData_printmyaddre(double* device_ptr);

}

#endif
