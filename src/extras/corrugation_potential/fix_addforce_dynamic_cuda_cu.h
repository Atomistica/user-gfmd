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

/*Copying from fix_viscous_cyda_cu.h Tristan Sharp*/

#include "cuda_shared.h"

extern "C" void Cuda_FixAddForce_DynamicCuda_Init(cuda_shared_data* sdata);
extern "C" void Cuda_FixAddForce_DynamicCuda_PostForce(cuda_shared_data* sdata, int groupbit,
	    	    F_FLOAT reflxcenter_sigma, F_FLOAT reflycenter_sigma,
	            F_FLOAT radspersigma_x, F_FLOAT radspersigma_y,
		    F_FLOAT aradsq, F_FLOAT border_dist,
		    F_FLOAT r_overEstar);
//extern "C" void Cuda_FixAddForceDynamicCuda_ComputeVector(cuda_shared_data* sdata);
