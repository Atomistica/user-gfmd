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
#include <stdlib.h>

#include "geometry.h"

/*!
 * Return a random number between 0 and 1.
 */
double rand1()
{
    return double(rand())/RAND_MAX;
}


/*!
 * Compute angle between two lines, counted clockwise.
 * I.e. if x1,y1 is rotated by the resulting angle it will lie on top of the
 * line to x2,y2.
 */
double clockwise_angle(double x1, double y1, double x2, double y2)
{
    double l1 = sqrt(x1*x1+y1*y1);
    double l2 = sqrt(x2*x2+y2*y2);
    
    /* dot-product */
    double cos_alpha = (x1*x2+y1*y2)/(l1*l2);
    /* cross-product */
    double sin_alpha = (x2*y1-x1*y2)/(l1*l2);
    
    return atan2(sin_alpha, cos_alpha);
}


/*!
 * Check if lines A1x,A1y-A2x,A2y and B1x,B1y-B2x,B2y cross
 */
bool lines_cross(double A1x, double A1y, double A2x, double A2y,
                 double B1x, double B1y, double B2x, double B2y)
{
    /* Cross product A1,A2 with B1,A2 */
    double c1 = (A1x-A2x)*(B1y-A2y)-(A1y-A2y)*(B1x-A2x);
    /* Cross product A1,A2 with B2,A2 */
    double c2 = (A1x-A2x)*(B2y-A2y)-(A1y-A2y)*(B2x-A2x);

    /* Both points B1 and B2 are on the same side of the (infinite) line that
       cuts through A1 and A2 */    
    if (c1*c2 > 0) {
        return false;
    }
    
    /* Cross product B1,B2 with A1,B2 */
    c1 = (B1x-B2x)*(A1y-B2y)-(B1y-B2y)*(A1x-B2x);
    /* Cross product B1,B2 with A2,B2 */
    c2 = (B1x-B2x)*(A2y-B2y)-(B1y-B2y)*(A2x-B2x);
    
    /* Both points A1 and A2 are on the same side of the (infinite) line that
       cuts through B1 and B2 */    
    if (c1*c2 > 0) {
        return false;
    }

    return true;
}


/*!
 * Check if point Bx,By is inside the triangle bounded by A1x,A1y, A2x,A2y
 * and A3x,A3y.
 */
bool tr_is_inside(double Bx, double By, double A1x, double A1y,
                  double A2x, double A2y, double A3x, double A3y)
{
    /* Cross product A1,A2 with B,A2 */
    double c1 = (A1x-A2x)*(By-A2y)-(A1y-A2y)*(Bx-A2x);
    /* Cross product A2,A3 with B,A3 */
    double c2 = (A2x-A3x)*(By-A3y)-(A2y-A3y)*(Bx-A3x);
    /* Cross product A3,A1 with B,A1 */
    double c3 = (A3x-A1x)*(By-A1y)-(A3y-A1y)*(Bx-A1x);

    if (c1*c2 > 0 && c2*c3 > 0 && c3*c1 > 0) {
        return true;
    }
    
    return false;
}


/*!
 * Compute area spanned by triangle A1x,A1y, A2x,A2y, A3x,A3y
 */
double tr_area(double A1x, double A1y, double A2x, double A2y, double A3x,
               double A3y)
{
    return 0.5*fabs((A1x-A2x)*(A3y-A2y)-(A1y-A2y)*(A3x-A2x));
}


/*!
 * Check if triangle bounded by Bx,By, A1x,A1y, A2x,A2y overlaps with the
 * triangle bounded by Cx,Cy, A1x,A1y, A2x,A2y.
 */
bool tr_overlap(double Bx, double By, double Cx, double Cy,
                double A1x, double A1y, double A2x, double A2y)
{
    /* Check if one triangle is fully inside the other one. A sufficient
       condition is that one point is inside the other triangle. */
    if (tr_is_inside(Bx,By,  Cx,Cy, A1x,A1y, A2x,A2y) ||
        tr_is_inside(Cx,Cy,  Bx,By, A1x,A1y, A2x,A2y)) {
        return true;
    }

    /* Check if triangle edges cross */
    if (lines_cross(A1x,A1y, Bx,By,  A2x,A2y, Cx,Cy) ||
        lines_cross(A2x,A2y, Bx,By,  A1x,A1y, Cx,Cy)) {
        return true;
    }

    /* Otherwise, triangles do not overlap */
    return false;
}


