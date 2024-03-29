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

#ifndef __PY_SURFACE_STIFFNESS_H
#define __PY_SURFACE_STIFFNESS_H

#include <Python.h>

#include "domain.h"
#include "error.h"
#include "memory.h"
#include "surface_stiffness.h"

typedef struct {
  PyObject_HEAD

  LAMMPS_NS::Error *error_;
  LAMMPS_NS::Memory *memory_;
  LAMMPS_NS::Domain *domain_;
  LAMMPS_NS::Force *force_;

  LAMMPS_NS::StiffnessKernel *kernel_;

} surface_stiffness_t;

extern PyTypeObject surface_stiffness_type;

#endif
