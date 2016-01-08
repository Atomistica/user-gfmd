
/*Copying from fix_viscous_cyda_cu.h Tristan Sharp*/

#include "cuda_shared.h"

extern "C" void Cuda_FixAddForce_DynamicCuda_Init(cuda_shared_data* sdata);
extern "C" void Cuda_FixAddForce_DynamicCuda_PostForce(cuda_shared_data* sdata, int groupbit,
	    	    F_FLOAT reflxcenter_sigma, F_FLOAT reflycenter_sigma,
	            F_FLOAT radspersigma_x, F_FLOAT radspersigma_y,
		    F_FLOAT aradsq, F_FLOAT border_dist,
		    F_FLOAT r_overEstar);
//extern "C" void Cuda_FixAddForceDynamicCuda_ComputeVector(cuda_shared_data* sdata);
