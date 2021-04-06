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
#include "atom.h"

using namespace LAMMPS_NS;

Atom::Atom(bigint in_natoms, Memory *in_memory)
{
  memory = in_memory;

  natoms = in_natoms;
  nlocal = natoms;
  nghost = 0;
  nmax = natoms;

  memory->create(x, nmax, 3, "Atom::x");
  memory->create(v, nmax, 3, "Atom::v");
  memory->create(f, nmax, 3, "Atom::f");
  memory->create(xeq, nmax, 3, "Atom::xeq");
  memory->create(mask, nmax, "Atom::mask");

  for (int i = 0; i < nmax; i++) {
    mask[i] = 1;
  }
}


Atom::~Atom()
{
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);
  memory->destroy(xeq);
}
