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
#ifndef __INTERPOLATED_MAP_H
#define __INTERPOLATED_MAP_H

#include <Python.h>

#include "error.h"
#include "memory.h"
#include "table2d.h"

typedef struct {
  PyObject_HEAD

  LAMMPS_NS::Error *error_;
  LAMMPS_NS::Memory *memory_;

  LAMMPS_NS::Table2D *map_;

} bicubic_t;

#endif
