/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

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
#include <math.h>
#include <string.h>

#include "table2d.h"

namespace LAMMPS_NS {

/*
 * values are supposed to be of size [0:nx][0:ny]
 */
Table2D::Table2D(int nx, int ny, double **values, bool interp, bool lowmem,
		 Error *error, Memory *memory)
{
  const int ix1[NCORN] = { 0,1,1,0 };
  const int ix2[NCORN] = { 0,0,1,1 };

  error_ = error;
  memory_ = memory;

  /*
   * calculate 2-d cubic parameters within each box.
   *
   * normalised coordinates.
   *  4--<--3
   *  |     ^
   *  v     |
   *  1-->--2
   */

  int irow, icol, nx1, nx2, npow1, npow2, npow1m, npow2m;
  int i, j;

  /* --- */

  interp_ = interp;

  nx_ = nx;
  ny_ = ny;
  int nboxs = nx_*ny_;

  /*
   * if lowmem = true then spline coefficients will be computed each
   * time eval is called
   */

  if (lowmem) {
    coeff_ = NULL;
    memory_->create(coeff_lowmem_, 4, 4, "Table2D::coeff_lowmem");
    values_ = values;
  }
  else {
    values_ = NULL;
    coeff_lowmem_ = NULL;
    memory_->create(coeff_, nboxs, 4, 4, "Table2D::coeff");
  }

  /*
   * for each box, create and solve the matrix equatoion.
   *    / values of  \     /              \     / function and \
   *  a |  products  | * x | coefficients | = b |  derivative  |
   *    \within cubic/     \ of 2d cubic  /     \    values    /
   */


  /*
   * construct the matrix.
   * this is the same for all boxes as coordinates are normalised.
   * loop through corners.
   */

  for (irow = 0; irow < NCORN; irow++) {
    nx1   = ix1[irow];
    nx2   = ix2[irow];
    /* loop through powers of variables. */
    for (npow1 = 0; npow1 <= 3; npow1++) {
      for (npow2 = 0; npow2 <= 3; npow2++) {
                         npow1m = npow1-1;
        if (npow1m < 0)  npow1m = 0;
                         npow2m = npow2-1;
        if (npow2m < 0)  npow2m=0;

        icol = 4*npow1+npow2;

        /* values of products within cubic and derivatives. */
        A_[irow   ][icol] = 1.0*(      pow(nx1,npow1 )      *pow(nx2,npow2 ) );
        A_[irow+4 ][icol] = 1.0*(npow1*pow(nx1,npow1m)      *pow(nx2,npow2 ) );
        A_[irow+8 ][icol] = 1.0*(      pow(nx1,npow1 )*npow2*pow(nx2,npow2m) );
        A_[irow+12][icol] = 1.0*(npow1*pow(nx1,npow1m)*npow2*pow(nx2,npow2m) );
      }
    }
  }


  /*
   * solve by gauss-jordan elimination with full pivoting.
   */

  GaussJordan(NPARA, A_[0], error_);

  if (coeff_) {

    for (int nhbox = 0; nhbox < nx_; nhbox++) {
      for (int ncbox = 0; ncbox < ny_; ncbox++) {
        int icol = ny_*nhbox+ncbox;
        compute_spline_coefficients(nhbox, ncbox, values, coeff_[icol]);
      }
    }

  }
}


Table2D::~Table2D()
{
  if (coeff_)
    memory_->destroy(coeff_);
  if (coeff_lowmem_)
    memory_->destroy(coeff_lowmem_);
}


void Table2D::eval(double x, double y, double &f)
{
  int r1box = (int) floor( x );
  int r2box = (int) floor( y );

  /*
   * find which box we're in and convert to normalised coordinates.
   */
  double x1   = x - r1box;
  double x2   = y - r2box;
  WRAPX(r1box);
  WRAPY(r2box);

  /*
   * get spline coefficients
   */
  double **coeffi = get_spline_coefficients(r1box, r2box);

  /*
   * compute splines
   */
  f    = 0.0;
  for (int i = 3; i >= 0; i--) {
    double sf = 0.0;
    for (int j = 3; j >= 0; j--) {
      sf = sf*x2 + coeffi[i][j];
    }
    f = f*x1 + sf;
  }
}


void Table2D::eval(double x, double y, double &f, double &dfdx, double &dfdy)
{
  int r1box = (int) floor( x );
  int r2box = (int) floor( y );

  /*
   * find which box we're in and convert to normalised coordinates.
   */
  double dx   = x - r1box;
  double dy   = y - r2box;
  WRAPX(r1box);
  WRAPY(r2box);

  /*
   * get spline coefficients
   */
  double **coeffi = get_spline_coefficients(r1box, r2box);

  /*
   * compute splines
   */
  f    = 0.0;
  dfdx = 0.0;
  dfdy = 0.0;
  for (int i = 3; i >= 0; i--) {
    double sf   = 0.0;
    double sfdy = 0.0;
    for (int j = 3; j >= 0; j--) {
      double      coefij = coeffi[i][j];
                  sf   =   sf*dy +   coefij;
      if (j > 0)  sfdy = sfdy*dy + j*coefij;
    }
                f    = f   *dx +   sf;
    if (i > 0)  dfdx = dfdx*dx + i*sf;
                dfdy = dfdy*dx +   sfdy;
  }
}


void Table2D::eval(double x, double y, double &f,
		   double &dfdx, double &dfdy,
		   double &d2fdxdx, double &d2fdydy, double &d2fdxdy)
  
{
  int r1box = (int) floor( x );
  int r2box = (int) floor( y );

  /*
   * find which box we're in and convert to normalised coordinates.
   */
  double dx   = x - r1box;
  double dy   = y - r2box;
  WRAPX(r1box);
  WRAPY(r2box);

  /*
   * get spline coefficients
   */
  double **coeffi = get_spline_coefficients(r1box, r2box);

  /*
   * compute splines
   */
  f       = 0.0;
  dfdx    = 0.0;
  dfdy    = 0.0;
  d2fdxdx = 0.0;
  d2fdydy = 0.0;
  d2fdxdy = 0.0;
  for (int i = 3; i >= 0; i--) {
    double sf      = 0.0;
    double sfdy    = 0.0;
    double s2fdydy = 0.0;
    double s2fdydx = 0.0;
    for (int j = 3; j >= 0; j--) {
      double      coefij  = coeffi[i][j];
                  sf      =   sf   *dy +         coefij;
      if (j > 0)  sfdy    = sfdy   *dy + j*      coefij;
      if (j > 1)  s2fdydy = s2fdydy*dy + j*(j-1)*coefij;
    }
                  f       = f      *dx +         sf;
    if (i > 0)    dfdx    = dfdx   *dx + i*      sf;
    if (i > 1)    d2fdxdx = d2fdxdx*dx + i*(i-1)*sf;
                  dfdy    = dfdy   *dx +         sfdy;
    if (i > 0)    d2fdxdy = d2fdxdy*dx + i*      sfdy;
  }
}

}
