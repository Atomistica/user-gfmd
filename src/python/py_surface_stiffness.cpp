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

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MATSCIPY_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <math.h>

#include "py_surface_stiffness.h"

#include "linearalgebra.h"

/* Allocate new instance */

PyObject *
surface_stiffness_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  surface_stiffness_t *self;

  self = (surface_stiffness_t *)type->tp_alloc(type, 0);

  self->error_ = NULL;
  self->memory_ = NULL;
  self->domain_ = NULL;
  self->force_ = NULL;
  self->kernel_ = NULL;

  return (PyObject *) self;
}


/* Release allocated memory */

void
surface_stiffness_dealloc(surface_stiffness_t *self)
{
  if (self->kernel_)
    delete self->kernel_;
  if (self->domain_)
    delete self->domain_;
  if (self->memory_)
    delete self->memory_;
  if (self->error_)
    delete self->error_;
  if (self->force_)
    delete self->force_;

  Py_TYPE(self)->tp_free((PyObject*) self);
}


/* Initialize instance */

#define MAX_ARG 1024

int
surface_stiffness_init(surface_stiffness_t *self, PyObject *args,
		      PyObject *kwargs)
{
  int narg, carg, i;
  double sx, sy;
  char *kernel_str, *kernel_args_in, *kernel_args;
  char *arg[MAX_ARG];

  if (!PyArg_ParseTuple(args, "ddss", &sx, &sy, &kernel_str, &kernel_args_in))
    return -1;

  kernel_args = strdup(kernel_args_in);

  narg = 0;
  i = 0;
  while (kernel_args[i] != '\0') {
    if (isspace(kernel_args[i])) {
      kernel_args[i] = '\0';
      i++;
    }
    while (isspace(kernel_args[i]))  i++;
    arg[narg] = &kernel_args[i];
    narg++;
    while (!isspace(kernel_args[i]) && kernel_args[i] != '\0')  i++;
  }
  
  self->error_ = new LAMMPS_NS::Error();
  self->memory_ = new LAMMPS_NS::Memory(self->error_);
  self->domain_ = new LAMMPS_NS::Domain();
  self->domain_->set_cell(sx, sy, 1.0);
  self->force_ = new LAMMPS_NS::Force();

  carg = 0;
  self->kernel_ = LAMMPS_NS::stiffness_kernel_factory(kernel_str, narg, &carg,
                                                      arg, self->domain_,
                                                      self->force_,
                                                      self->memory_,
                                                      self->error_);

  free(kernel_args);

  return 0;
}


/* Call object */

PyObject *
surface_stiffness_call(surface_stiffness_t *self, PyObject *args,
		      PyObject *kwargs)
{
  PyObject *py_qx, *py_qy, *py_dU = NULL;

  /* We support passing coordinates (x, y), numpy arrays (x, y)
     and numpy arrays r */

  py_qx = NULL;
  py_qy = NULL;
  if (!PyArg_ParseTuple(args, "OO|O", &py_qx, &py_qy, &py_dU))
    return NULL;

  /* This should a single numpy array r */

  py_qx = PyArray_FROMANY(py_qx, NPY_DOUBLE, 0, 2, NPY_C_CONTIGUOUS);
  if (!py_qx)
    return NULL;
  py_qy = PyArray_FROMANY(py_qy, NPY_DOUBLE, 0, 2, NPY_C_CONTIGUOUS);
  if (!py_qy)
    return NULL;
  if (py_dU) {
    py_dU = PyArray_FROMANY(py_dU, NPY_COMPLEX128, 0, 2, NPY_C_CONTIGUOUS);
    if (!py_dU)
      return NULL;
  }

  npy_intp ndims = PyArray_NDIM(py_qx);
  if (PyArray_NDIM(py_qy) != ndims) {
    PyErr_SetString(PyExc_RuntimeError, "qx and qy arrays need to have "
		    "same shape.");
    return NULL;
  }
  if (py_dU) {
    if (PyArray_NDIM(py_dU) != ndims) {
      PyErr_SetString(PyExc_RuntimeError, "qx, qy and dU arrays need to have "
		      "same shape.");
      return NULL;
    }
  }

  npy_intp *dims = PyArray_DIMS(py_qx);
  npy_intp exdims[16];
  int n = 1;
  for (int i = 0; i < ndims; i++) {
    n *= dims[i];
    exdims[i] = dims[i];

    if (PyArray_DIM(py_qy, i) != dims[i]) {
      PyErr_SetString(PyExc_RuntimeError, "qx and qy arrays need to have "
		      "same shape.");
      return NULL;
    }
    if (py_dU) {
      if (PyArray_DIM(py_dU, i) != dims[i]) {
	PyErr_SetString(PyExc_RuntimeError, "qx, qy and dU arrays need to have "
			"same shape.");
	return NULL;
      }
    }
  }

  exdims[ndims] = self->kernel_->get_dimension();
  exdims[ndims+1] = self->kernel_->get_dimension();

  double *qx = (double *) PyArray_DATA(py_qx);
  double *qy = (double *) PyArray_DATA(py_qy);
  double_complex *dU = NULL;
  if (py_dU) {
    dU = (double_complex *) PyArray_DATA(py_dU);
  }

  int phi_off = exdims[ndims]*exdims[ndims+1];

  PyObject *py_phi = PyArray_SimpleNew(ndims+2, exdims, NPY_COMPLEX128);
  double_complex *phi = (double_complex *) PyArray_DATA(py_phi);

  for (int i = 0; i < n; i++) {
    if (dU) {
      self->kernel_->get_stiffness_matrix(*qx, *qy, phi, *dU);
      dU++;
    }
    else {
      self->kernel_->get_stiffness_matrix(*qx, *qy, phi);
    }
    qx++;
    qy++;
    phi += phi_off;
  }

  Py_DECREF(py_qx);
  Py_DECREF(py_qy);
  if (py_dU) {
    Py_DECREF(py_dU);
  }

  return py_phi;
}


/* Set model parameters */

PyObject *
surface_stiffness_set_parameters(surface_stiffness_t *self, PyObject *args,
				 PyObject *kwargs)
{
  Py_ssize_t size = PyTuple_GET_SIZE(args);

  if (size != self->kernel_->get_number_of_parameters()) {
    char errstr[1024];
    sprintf(errstr, "%i arguments provided in call to set_parameters, but "
	    "kernel takes %i parameters.", (int) size,
	    self->kernel_->get_number_of_parameters());
    PyErr_SetString(PyExc_RuntimeError, errstr);
    return NULL;
  }

  double pars[size];

  for (Py_ssize_t i = 0; i < size; i++) {
    pars[i] = PyFloat_AsDouble(PyTuple_GET_ITEM(args, i));
    if (PyErr_Occurred()) {
      return NULL;
    }
  }

  self->kernel_->set_parameters(pars);

  Py_RETURN_NONE;
}


/* Methods declaration */

PyMethodDef surface_stiffness_methods[] = {
  { "set_parameters",
    (PyCFunction) surface_stiffness_set_parameters, METH_VARARGS,
    "Set model parameters." },
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


/* Class declaration */

PyTypeObject surface_stiffness_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "surface_stiffness.StiffnessKernel",        /* tp_name */
    sizeof(surface_stiffness_t),                /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)surface_stiffness_dealloc,      /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    (ternaryfunc)surface_stiffness_call,        /* tp_call */
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
    surface_stiffness_methods,                  /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)surface_stiffness_init,           /* tp_init */
    0,                                          /* tp_alloc */
    surface_stiffness_new,                      /* tp_new */
};

