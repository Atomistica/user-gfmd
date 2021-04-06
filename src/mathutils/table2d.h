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
#ifndef TABLE2D_H
#define TABLE2D_H

#include "error.h"
#include "memory.h"

#include "linearalgebra.h"

#define NPARA (4*4)   // 4^dim
#define NCORN 4

#define WRAPX(x) { while (x >= nx_) x -= nx_; while (x < 0) x += nx_; }
#define WRAPY(y) { while (y >= ny_) y -= ny_; while (y < 0) y += ny_; }

namespace LAMMPS_NS {

class Table2D {
 public:
  Table2D(int, int, double **, bool, bool, Error *, Memory *);
  ~Table2D();

  void eval(double, double, double &);
  void eval(double, double, double &, double &, double &);
  void eval(double, double, double &, double &, double &, double &, double &,
            double &);

 protected:
  /* LAMMPS stuff */
  Error *error_;
  Memory *memory_;

  /* table dimensions */
  int nx_, ny_;

  /* interpolate derivatives */
  bool interp_;

  /* values */
  double **values_;

  /* spline coefficients */
  double ***coeff_;

  /* spline coefficients if lowmem is true */
  double **coeff_lowmem_;

  /* lhs matrix */
  double A_[NPARA][NPARA];

  double **get_spline_coefficients(int r1box, int r2box) {
    if (coeff_) {
      int ibox = ny_*r1box+r2box;
      return coeff_[ibox];
    }
    else {
      compute_spline_coefficients(r1box, r2box, values_, coeff_lowmem_);
      return coeff_lowmem_;
    }
  }

  void compute_spline_coefficients(int nhbox, int ncbox, double **values,
				   double **coeff) {
    const int ix1[NCORN] = { 0,1,1,0 };
    const int ix2[NCORN] = { 0,0,1,1 };

    /*
     * construct the 16 r.h.s. vectors ( 1 for each box ).
     * loop through boxes.
     */ 

    double B[NPARA];

    for (int irow = 0; irow < NCORN; irow++) {
      int nx1  = ix1[irow]+nhbox;
      int nx2  = ix2[irow]+ncbox;
      /* wrap to box */
      WRAPX(nx1);
      WRAPY(nx2);
      /* values of function and derivatives at corner. */
      B[irow   ] = values[nx1][nx2];
      /* interpolate derivatives */
      if (interp_) {
        int nx1p = nx1+1;
        int nx1m = nx1-1;
        int nx2p = nx2+1;
        int nx2m = nx2-1;
        WRAPX(nx1p);
        WRAPX(nx1m);
        WRAPY(nx2p);
        WRAPY(nx2m);
        B[irow+4 ] = (values[nx1p][nx2]-values[nx1m][nx2])/2;
        B[irow+8 ] = (values[nx1][nx2p]-values[nx1][nx2m])/2;
      }
      else {
        B[irow+4 ] = 0.0;
        B[irow+8 ] = 0.0;
      }
      B[irow+12] = 0.0;
    }

    double tmp[NPARA];
    mat_mul_vec(NPARA, &A_[0][0], B, tmp);

    /*
     * get the coefficient values.
     */

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        int irow = 4*i+j;
	
        coeff[i][j] = tmp[irow];
      }
    }
  }
};

}

#endif
