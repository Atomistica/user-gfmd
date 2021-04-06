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
#ifndef ATOM_VEC_GFMD_CUDA_CU_H_
#define ATOM_VEC_GFMD_CUDA_CU_H_

extern "C" {

void Cuda_AtomVecGFMDCuda_Init(cuda_shared_data* sdata);
int Cuda_AtomVecGFMDCuda_PackExchangeList(cuda_shared_data* sdata, int n,
					  int dim, void* buf_send);
int Cuda_AtomVecGFMDCuda_PackExchange(cuda_shared_data* sdata, int nsend,
				      void* buf_send, void* copylist);
int Cuda_AtomVecGFMDCuda_UnpackExchange(cuda_shared_data* sdata, int nsend,
					void* buf_send, void* copylist);
int Cuda_AtomVecGFMDCuda_PackBorder(cuda_shared_data* sdata, int nsend,
				    int iswap, void* buf_send, int* pbc, 
				    int pbc_flag);
int Cuda_AtomVecGFMDCuda_PackBorderVel(cuda_shared_data* sdata, int n,
				       int iswap, void* buf_send, int* pbc, 
				       int pbc_flag);
int Cuda_AtomVecGFMDCuda_PackBorder_Self(cuda_shared_data* sdata, int n,
					 int iswap, int first, int* pbc,
					 int pbc_flag);
int Cuda_AtomVecGFMDCuda_PackBorderVel_Self(cuda_shared_data* sdata, int n,
					    int iswap, int first, int* pbc,
					    int pbc_flag);
int Cuda_AtomVecGFMDCuda_UnpackBorder(cuda_shared_data* sdata, int n, int first,
				      void* buf_recv);
int Cuda_AtomVecGFMDCuda_UnpackBorderVel(cuda_shared_data* sdata, int n,
					 int first, void* buf_recv);

}

#endif /*ATOM_VEC_GFMD_CUDA_CU_H_*/
