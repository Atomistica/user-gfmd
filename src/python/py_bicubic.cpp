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
/*
 * A wrapper around Table2D to be used from within Python+Numpy
 */

#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>

#include "py_bicubic.h"

/* Allocate new instance */

static PyObject *
bicubic_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  bicubic_t *self;

  self = (bicubic_t *)type->tp_alloc(type, 0);

  self->error_ = NULL;
  self->memory_ = NULL;
  self->map_ = NULL;

  return (PyObject *) self;
}


/* Release allocated memory */

static void
bicubic_dealloc(bicubic_t *self)
{
  if (self->map_)
    delete self->map_;
  if (self->memory_)
    delete self->memory_;
  if (self->error_)
    delete self->error_;

  self->ob_type->tp_free((PyObject*) self);
}


/* Initialize instance */

static int
bicubic_init(bicubic_t *self, PyObject *args,
		      PyObject *kwargs)
{
  PyObject *py_map_data_in;

  if (!PyArg_ParseTuple(args, "O|O", &py_map_data_in))
    return -1;

  PyObject *py_map_data;
  npy_intp nx, ny;

  py_map_data = PyArray_FROMANY(py_map_data_in, NPY_DOUBLE, 2, 2, 0);
  if (!py_map_data)
    return -1;
  nx = PyArray_DIM(py_map_data, 0);
  ny = PyArray_DIM(py_map_data, 1);

  self->error_ = new LAMMPS_NS::Error();
  self->memory_ = new LAMMPS_NS::Memory(self->error_);

  double **map_data = new double*[nx];
  for (int i = 0; i < nx; i++) {
    map_data[i] = &((double *) PyArray_DATA(py_map_data))[i*ny];
  }

  self->map_ = new LAMMPS_NS::Table2D(nx, ny, map_data, true, false,
                                      self->error_, self->memory_);

  delete map_data;
  Py_DECREF(py_map_data);

  return 0;
}


/* Call object */

static PyObject *
bicubic_call(bicubic_t *self, PyObject *args,
		      PyObject *kwargs)
{
  PyObject *py_x, *py_y;

  /* We support passing coordinates (x, y), numpy arrays (x, y)
     and numpy arrays r */

  py_x = NULL;
  py_y = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &py_x, &py_y))
    return NULL;

  if (!py_y) {
    /* This should a single numpy array r */

    PyObject *py_r;
    py_r = PyArray_FROMANY(py_x, NPY_DOUBLE, 2, 2, 0);
    if (!py_r)
      return NULL;

    if (PyArray_DIM(py_r, 1) != 2) {
      PyErr_SetString(PyExc_TypeError, "Map index needs to have x- and y- "
		      "component only.");
      return NULL;
    }

    npy_intp n = PyArray_DIM(py_r, 0);
    double *r = (double *) PyArray_DATA(py_r);

    PyObject *py_v = PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    double *v = (double *) PyArray_DATA(py_v);

    for (int i = 0; i < n; i++) {
      double dx, dy;
      self->map_->eval(r[2*i], r[2*i+1], v[i], dx, dy);
    }

    Py_DECREF(py_r);

    return py_v;
  }
  else if (( PyFloat_Check(py_x)||PyInt_Check(py_x)||PyLong_Check(py_x) ) &&
	   ( PyFloat_Check(py_y)||PyInt_Check(py_y)||PyLong_Check(py_y) )) {
    /* x and y are specified separately, and are scalars */

    double v, dx, dy;
    self->map_->eval(PyFloat_AsDouble(py_x), PyFloat_AsDouble(py_y),
                     v, dx, dy);
    return PyFloat_FromDouble(v);
  }
  else {
    /* x and y are specified separately */
    PyObject *py_xd, *py_yd;
    py_xd = PyArray_FROMANY(py_x, NPY_DOUBLE, 1, 3, 0);
    if (!py_xd)
      return NULL;
    py_yd = PyArray_FROMANY(py_y, NPY_DOUBLE, 1, 3, 0);
    if (!py_yd)
      return NULL;

    /* Check that x and y have the same number of dimensions */
    if (PyArray_NDIM(py_xd) != PyArray_NDIM(py_yd)) {
      PyErr_SetString(PyExc_TypeError, "x- and y-components need to have "
		      "identical number of dimensions.");
      return NULL;
    }

    /* Check that x and y have the same length in each dimension */
    int ndims = PyArray_NDIM(py_xd);
    npy_intp *dims = PyArray_DIMS(py_xd);
    npy_intp n = 1;
    for (int i = 0; i < ndims; i++) {
      npy_intp d = PyArray_DIM(py_yd, i);

      if (dims[i] != d) {
	PyErr_SetString(PyExc_TypeError, "x- and y-components vectors need to "
			"have the same length.");
	return NULL;
      }

      n *= d;
    }

    double *x = (double *) PyArray_DATA(py_xd);
    double *y = (double *) PyArray_DATA(py_yd);

    PyObject *py_v = PyArray_SimpleNew(ndims, dims, NPY_DOUBLE);
    double *v = (double *) PyArray_DATA(py_v);

    for (int i = 0; i < n; i++) {
      double dx, dy;
      self->map_->eval(x[i], y[i], v[i], dx, dy);
    }

    Py_DECREF(py_xd);
    Py_DECREF(py_yd);

    return py_v;
  }
}


/* Class declaration */

PyTypeObject bicubic_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /* ob_size */
    "bicubic.Bicubic",                          /* tp_name */
    sizeof(bicubic_t),                          /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)bicubic_dealloc,                /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    (ternaryfunc)bicubic_call,                  /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Potential objects",                        /* tp_doc */
    0,		                                /* tp_traverse */
    0,		                                /* tp_clear */
    0,		                                /* tp_richcompare */
    0,		                                /* tp_weaklistoffset */
    0,		                                /* tp_iter */
    0,		                                /* tp_iternext */
    0,                                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)bicubic_init,                     /* tp_init */
    0,                                          /* tp_alloc */
    bicubic_new,                                /* tp_new */
};


extern "C"
void initbicubic()  {
  PyObject *m;

  m = Py_InitModule("bicubic", NULL);
  if (!m)
    return;

  import_array();

  if (PyType_Ready(&bicubic_type) < 0)
    return;

  Py_INCREF(&bicubic_type);
  PyModule_AddObject(m, "Bicubic",
		     (PyObject *) &bicubic_type);
}
