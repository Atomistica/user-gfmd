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
