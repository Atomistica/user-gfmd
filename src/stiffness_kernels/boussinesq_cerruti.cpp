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
/*
 * Boussinesq-Cerruti solution for uniform pressure on a triangular area.
 * Formulae and notation is from: Li, Berger, J. Elast. 63, 137 (2001)
 *
 * Sign convention: negative P is pushing down, positive P is pushing up.
 * Similarly: negative uz is displacement into the surface.
 * I.e. signs are reverse to what is usually assumed (Johnson; Li, Berger)
 */

#include <math.h>

#include "boussinesq_cerruti.h"

/*!
 * Return displacements ux,uy,uz at position x,y,z due to a force P,Qx,Qy at
 * the origin.
 * Note: negativ z is within the solid, and negative P is pushing down.
 */
void pt_PQ(double G, double nu,                  /* elastic moduli */
           double P, double Qx, double Qy,       /* force */
           double x, double y, double z,         /* position */
           double &ux, double &uy, double &uz)   /* displacement at position */
{
    double fac = 1/(4*M_PI*G);
    double a = 1-2*nu;
    double b = 2*(1-nu);
    
    double r_sq = x*x+y*y;
    double r = sqrt(r_sq);
    double r_3 = r*r_sq;
    
    double rz = r-z;
    double rz_sq = rz*rz;
    
    /*
     * Normal force - Boussinesq solution
     * See Johnson p.50
     */
    ux += -fac*P*( -x*z/r_3 - a*x/(r*rz) );
    uy += -fac*P*( -y*z/r_3 - a*y/(r*rz) );
    uz +=  fac*P*(  z*z/r_3 + b/r );
    
    /*
     * Tangential tractions - Cerruti solution
     * See Johnson p.69
     */
    ux +=  fac*Qx*(  1.0/r + x*x/(r*r_sq)+a*( 1/rz-x*x/(r*r_sq) ) );
    uy +=  fac*Qx*(  x*y/r_3 - a*x*y/(r*rz_sq) );
    uz += -fac*Qx*( -x*z/r_3 + a*x/(r*rz) );

    ux +=  fac*Qy*(  x*y/r_3 - a*x*y/(r*rz_sq) );
    uy +=  fac*Qy*(  1.0/r + y*y/(r*r_sq)+a*( 1/rz-y*y/(r*r_sq) ) );
    uz += -fac*Qy*( -y*z/r_3 + a*y/(r*rz) );
}
