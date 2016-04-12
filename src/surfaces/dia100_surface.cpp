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

#include "dia100_surface.h"

/*!
 * Constructor that accepts LAMMPS string arguments
 */
Diamond100Surface::Diamond100Surface(int narg, int *carg, char **arg,
                                     Error *error)
    : CrystalSurface()
{
    strcpy(name_, "dia100");

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
Diamond100Surface::Diamond100Surface(double lattice_constant, int nu,
                                     Error *error)
    : CrystalSurface(), lattice_constant_(lattice_constant), nu_(nu)
{
    strcpy(name_, "dia100");

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
void Diamond100Surface::dump_info(FILE *f)
{
    CrystalSurface::dump_info(f);

    fprintf(f, "Lattice constant = %f; primitive unit cell is repeated %i "
            "times in the surface normal direction.\n", lattice_constant_, nu_);
}


/*!
 * Create primitive unit cell
 */
void Diamond100Surface::primitive_cell(Error *error)
{
    // a0 is the conventional cubic lattice constant (in units of nearest
    // neighbor distance)
    double a0 = 4/sqrt(3);

    // cell vectors (in units of nearest neighbor distances)
    // v    v    v - vectors are columns here
    double _cell[9] = {
      sqrt(2)*a0/2,             0,  sqrt(2)*a0/4,
                 0,  sqrt(2)*a0/2,  sqrt(2)*a0/4,
                 0,             0,          a0/2
    };
    set_cell(_cell, error);

    // two atoms per surface unit cell
    num_atoms_ = 2;

    // atom positions
    double positions[6] = {
                 0,  0,     0,
      sqrt(2)*a0/4,  0,  a0/4
    };
    positions_ = new double[6];
    std::copy(positions, positions+6, positions_);

    // all types initially 1
    types_ = new int[2];
    types_[0] = 1;
    types_[1] = 1;
}
