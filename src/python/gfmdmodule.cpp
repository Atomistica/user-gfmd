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

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MATSCIPY_ARRAY_API
#include <numpy/arrayobject.h>

#include <stdbool.h>
#include <stddef.h>

#include "py_bicubic.h"
#include "py_surface_stiffness.h"

#include "gfmdmodule.h"

/*
 * Method declaration
 */

static PyMethodDef module_methods[] = {
    { NULL, NULL, 0, NULL }  /* Sentinel */
};

/*
 * Module initialization
 */

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

/*
 * Module declaration
 */

#if PY_MAJOR_VERSION >= 3
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
        ob = PyModule_Create(&moduledef);
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        ob = Py_InitModule3(name, methods, doc);
#endif

extern "C"
MOD_INIT(_gfmd)
{
    PyObject* m;

    import_array();

    MOD_DEF(m, "_gfmd", module_methods,
            "Python interface to GFMD.");

    if (PyType_Ready(&bicubic_type) < 0)
       goto fail;

    Py_INCREF(&bicubic_type);
    PyModule_AddObject(m, "Bicubic",
                       (PyObject *) &bicubic_type);

    if (PyType_Ready(&surface_stiffness_type) < 0)
        goto fail;

    Py_INCREF(&surface_stiffness_type);
    PyModule_AddObject(m, "SurfaceStiffness",
                       (PyObject *) &surface_stiffness_type);

fail:
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
