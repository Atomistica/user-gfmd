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
#ifndef LMP_ATOM_H
#define LMP_ATOM_H

#include "lmptype.h"
#include "memory.h"

namespace LAMMPS_NS {

class Atom {
 public:
  bigint natoms;                // total # of atoms in system, could be 0
  int nlocal,nghost;            // # of owned and ghost atoms on this proc
  int nmax;                     // max # of owned+ghost in arrays on this proc

  // per-atom arrays
  // customize by adding new array

  int *tag,*type,*mask;
  double **x,**v,**f,**xeq;

  Atom(bigint, Memory *);
  virtual ~Atom();

 protected:
  Memory *memory;
};

}

#endif
