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
#include <math.h>

#include "linearalgebra.h"

#include "fcc111_surface.h"

/*!
 * Constructor, pre-compute neighbor list (currently up to 3rd neighbor shell)
 */
FCC111Surface::FCC111Surface(int narg, int *carg, char **arg, Error *error)
    : CrystalSurface()
{
    strcpy(name_, "fcc111");

    read_double(error, narg, arg, carg, &lattice_constant_, "lattice constant", 0, 100);

    read_integer(error, narg, arg, carg, &nu_, "number of lattice cells", 1, 100);

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
void FCC111Surface::dump_info(FILE *f)
{
    CrystalSurface::dump_info(f);

    fprintf(f, "Lattice constant = %f; primitive unit cell is repeated %i "
            "times in the surface normal direction.\n", lattice_constant_, nu_);
}


/*!
 * Create primitive unit cell
 */
void FCC111Surface::primitive_cell(Error *error)
{
    // Cell vectors (in units of nearest neighbor distances).
    //                  v  v  v  - vectors are columns here
    double _cell[9] = { 1, 0,         0.5,
                        0, sqrt(3.0), 1.0/(2.0*sqrt(3.0)),
                        0, 0,         sqrt(2.0/3.0) };
    set_cell(_cell, error);

    // Two atoms per surface unit cell.
    num_atoms_ = 2;

    // Atoms are located in the bottom left corner and in the middle of cell.
    double scaled_positions[6] = { 0,0,0, 0.5,0.5,0 };
    positions_ = new double[6];
    mat<double> cell(3, cell_);
    // Convert from scaled to unscaled positions.
    for (int i = 0; i < 2; i++) {
        vec<double> sr(3, &scaled_positions[3*i]);
        vec<double> r(3, &positions_[3*i]);
        r = cell.dot(sr);
    }

    // all types initially 1
    types_ = new int[2];
    types_[0] = 1;
    types_[1] = 1;
}
