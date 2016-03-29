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
