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
