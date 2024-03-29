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
/*
 * Boussinesq-Cerruti solution for uniform pressure on a square area.
 * Formulae and notation are from Pohrt, Qiang (Physical Mesomechanics,
 * to be published).
 *
 * Sign convention: negative P is pushing down, positive P is pushing up.
 * Similarly: negative uz is displacement into the surface.
 * I.e. signs are reverse to what is usually assumed (Johnson; Li, Berger)
 */

#include <assert.h> 
#include <math.h>
#include <stdio.h>

#include "geometry.h"

#include "pohrt_li.h"

namespace PohrtLi {

/*!
 * Compute displacements at point Bx, By due to uniform pressure P
 * and tractions Qx, Qy on a rectangular area centered at the origin with
 * sides of length 2a,2b.
 */
void sq_uniform_PQ(double G, double nu,
                   double P, double Qx, double Qy, double a, double b, 
                   double Bx, double By, double &ux, double &uy, double &uz)
{
    double k = Bx+a;
    double l = Bx-a;
    double m = By+b;
    double n = By-b;

    /*
     * Displacements due to normal load P
     */
    ux += (1-2*nu)*P/(4*M_PI*G)*
        ( 0.5*(m*log((k*k+m*m)/(l*l+m*m))+n*log((l*l+n*n)/(k*k+n*n))) +
        k*(atan(m/k)-atan(n/k)) + l*(atan(n/l)-atan(m/l)) );
    uy += (1-2*nu)*P/(4*M_PI*G)*
        ( 0.5*(k*log((k*k+m*m)/(k*k+n*n))+l*log((l*l+n*n)/(l*l+m*m))) +
        m*(atan(k/m)-atan(l/m)) + n*(atan(l/n)-atan(k/n)) );
    uz += (1-nu)*P/(2*M_PI*G)*
        ( k*log((m+sqrt(k*k+m*m))/(n+sqrt(k*k+n*n))) +
        l*log((n+sqrt(l*l+n*n))/(m+sqrt(l*l+m*m))) +
        m*log((k+sqrt(k*k+m*m))/(l+sqrt(l*l+m*m))) +
        n*log((l+sqrt(l*l+n*n))/(k+sqrt(k*k+n*n))) );

    /*
     * Displacements due to traction Qx
     */
    ux += 1.0*Qx/(2*M_PI*G)*
        ( (1-nu)*( k*log((m+sqrt(k*k+m*m))/(n+sqrt(k*k+n*n))) +
        l*log((n+sqrt(l*l+n*n))/(m+sqrt(l*l+m*m))) ) +
        m*log((k+sqrt(k*k+m*m))/(l+sqrt(l*l+m*m))) +
        n*log((l+sqrt(l*l+n*n))/(k+sqrt(k*k+n*n))) );
    uy += nu*Qx/(2*M_PI*G)*
        ( sqrt(n*n+k*k) - sqrt(m*m+k*k) +
        sqrt(m*m+l*l) - sqrt(n*n+l*l) );
    uz += -(1-2*nu)*Qx/(4*M_PI*G)*
        ( 0.5*(m*log((k*k+m*m)/(l*l+m*m))+n*log((l*l+n*n)/(k*k+n*n))) +
        k*(atan(m/k)-atan(n/k)) + l*(atan(n/l)-atan(m/l)) );

    /*
     * Displacements due to traction Qy
     */
    ux += nu*Qy/(2*M_PI*G)*
        ( sqrt(n*n+k*k) - sqrt(m*m+k*k) +
        sqrt(m*m+l*l) - sqrt(n*n+l*l) );
    uy += 1.0*Qy/(2*M_PI*G)*
        ( (1-nu)*( m*log((k+sqrt(m*m+k*k))/(l+sqrt(m*m+l*l))) +
        n*log((l+sqrt(n*n+l*l))/(k+sqrt(n*n+k*k))) ) +
        k*log((m+sqrt(m*m+k*k))/(n+sqrt(n*n+k*k))) +
        l*log((n+sqrt(n*n+l*l))/(m+sqrt(m*m+l*l))) );
    uz += -(1-2*nu)*Qy/(4*M_PI*G)*
        ( 0.5*(k*log((k*k+m*m)/(k*k+n*n))+l*log((l*l+n*n)/(l*l+m*m))) +
        m*(atan(k/m)-atan(l/m)) + n*(atan(l/n)-atan(k/n)) );
}

}
