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

#include "linearalgebra.h"

#include "dia111_surface.h"

/*!
 * Constructor that accepts LAMMPS string arguments
 */
Diamond111Surface::Diamond111Surface(int narg, int *carg, char **arg,
                                     Error *error)
    : CrystalSurface()
{
    strcpy(name_, "dia111");

    read_double(error, narg, arg, carg, &lattice_constant_, "lattice constant",
                0, 100);

    read_integer(error, narg, arg, carg, &nu_, "number of lattice cells", 1,
                 100);

    primitive_cell(error);
    if (nu_ > 1) {
        supercell(1, 1, nu_, error);
        reverse_indices();
    }

    mat<double> cell(3, cell_);
    set_cell_and_scale_atoms((lattice_constant_*cell).const_data(), error);
}


/*!
 * Constructor
 */
Diamond111Surface::Diamond111Surface(double lattice_constant, int nu,
                                     Error *error)
    : CrystalSurface(), lattice_constant_(lattice_constant), nu_(nu)
{
    strcpy(name_, "dia111");

    primitive_cell(error);
    if (nu_ > 1) {
        supercell(1, 1, nu_, error);
        reverse_indices();
    }

    mat<double> cell(3, cell_);
    set_cell_and_scale_atoms((lattice_constant_*cell).const_data(), error);
}


/*!
 * Dump generic text info to file
 */
void Diamond111Surface::dump_info(FILE *f)
{
    CrystalSurface::dump_info(f);

    fprintf(f, "Lattice constant = %f; primitive unit cell is repeated %i "
            "times in the surface normal direction.\n", lattice_constant_, nu_);
}


/*!
 * Create primitive unit cell
 */
void Diamond111Surface::primitive_cell(Error *error)
{
    // cell vectors (in units of nearest neighbor distances)
    double a1 = 2*sqrt(2);
    double a2 = 2*sqrt(2.0/3);
    double a3 = 4./3;
    // v    v    v - vectors are columns here
    double _cell[9] = {
        a1,   0,  1./3*a1,
         0,  a2,        0,
         0,   0,       a3
    };
    set_cell(_cell, error);

    // two atoms per surface unit cell
    num_atoms_ = 4;

    // atom positions
    double positions[12] = {
               0,       0,        0,
         1./3*a1,       0,  1./4*a3,
         1./2*a1,  0.5*a2,        0,
         5./6*a1,  0.5*a2,  1./4*a3
    };
    positions_ = new double[12];
    std::copy(positions, positions+12, positions_);

    // all types initially 1
    types_ = new int[4];
    std::fill(types_, types_+4, 1);
}
