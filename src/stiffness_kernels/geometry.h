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
#ifndef __GEOMETRY_H
#define __GEOMETRY_H

double rand1();

double clockwise_angle(double x1, double y1, double x2, double y2);

bool lines_cross(double A1x, double A1y, double A2x, double A2y,
                 double B1x, double B1y, double B2x, double B2y);

bool tr_is_inside(double Bx, double By, double A1x, double A1y,
                  double A2x, double A2y, double A3x, double A3y);
                  
double tr_area(double A1x, double A1y, double A2x, double A2y, double A3x,
               double A3y);

bool tr_overlap(double Bx, double By, double Cx, double Cy,
                double A1x, double A1y, double A2x, double A2y);

#endif
