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
#ifndef __LI_BERGER_H
#define __LI_BERGER_H

namespace LiBerger {

void tr_uniform_PQ_c(double G, double nu,
                     double P, double Qx, double Qy, double Bx, double By,
		     double A2x, double A2y, double A3x, double A3y,
		     double &ux, double &uy, double &uz);

void tr_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double A1x, double A1y,
		   double A2x, double A2y, double A3x, double A3y,
		   double Bx, double By, double &ux, double &uy, double &uz);
		           
void tr_linear_PQ_c(double G, double nu,
                    double P, double Qx, double Qy, double Bx, double By,
                    double A2x, double A2y, double A3x, double A3y,
                    double &ux, double &uy, double &uz);

void tr_linear_PQ(double G, double nu,
                  double P, double Qx, double Qy, double A1x, double A1y,
                  double A2x, double A2y, double A3x, double A3y,
                  double Bx, double By, double &ux, double &uy, double &uz);

void sq_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double a, double b, 
		   double Bx, double By, double &ux, double &uy, double &uz);

}

#endif
