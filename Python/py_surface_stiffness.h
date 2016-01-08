#ifndef __INTERPOLATED_MAP_H
#define __INTERPOLATED_MAP_H

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

  LAMMPS_NS::StiffnessKernel *kernel_;

} surface_stiffness_t;

#endif
