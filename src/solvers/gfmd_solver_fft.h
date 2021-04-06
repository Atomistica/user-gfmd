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
#ifndef GFMD_SOLVER_FFT_H
#define GFMD_SOLVER_FFT_H

#include "gfmd_solver.h"
#include "linearalgebra.h"
#include "pointers.h"

namespace LAMMPS_NS {

class GFMDSolverFFT : public GFMDSolver {
 public:
  GFMDSolverFFT(LAMMPS *);
  virtual ~GFMDSolverFFT();

  virtual void set_grid_size(int, int, int);

  virtual double memory_usage();

  virtual void dump_stiffness();
  virtual void dump_greens_function();

  virtual void fft_forward(double **, double_complex **, 
			   double_complex **fac=NULL);
  virtual void fft_reverse(double_complex **, double **,
			   double_complex **fac=NULL);

  double_complex **phi;

 protected:
  class FFT3d *fft;
  double *fft_data;

  void dump(char *, double_complex **);
};

}

#endif
