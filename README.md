Green's function molecular dynamics (GFMD) code
===============================================

This package implements the Green's function molecular dynamics, by introducing
a LAMMPS fix, `fix gfmd`. It is used to evaluate the effective forces acting on
an elastic half-space from the dynamical matrix.

This is an implementation of the method described in:

- [C. Campana, M.H. Muser, Phys. Rev. B 74, 075420 (2006)](http://dx.doi.org/10.1103/PhysRevB.74.075420)
- [L. Pastewka, T.A. Sharp, M.O. Robbins, Phys. Rev. B 86, 075459 (2012)](http://dx.doi.org/10.1103/PhysRevB.86.075459)

If you use this code, please cite these publications.

Build status
------------

[![Build Status](https://travis-ci.org/Atomistica/user-gfmd.svg?branch=master)](https://travis-ci.org/Atomistica/user-gfmd)

LAMMPS
------

The files in this directory contain a user-contributed package for LAMMPS. See
the documentation files of this commands for detailed usage. These can be
found in the `doc` subdirectory.

Compiling this package along with LAMMPS requires the FFT3d wrappers from the
kspace package of LAMMPS to be included.

You should also append 'user-gfmd' to the 'PACKUSER' variable in the file of
src/Makefile before invoke `make yes-user`. Note that `make yes-user-gfmd` does
not work if used with the USER-CUDA package. If used with USER-CUDA, go into
the USER-GFMD directory and execute `sh Install.sh 1` manually.

Functionality beyond pair potentials will be included soon.

Testing
-------

There are two separate testing frameworks:

1. Unittests are separate from the LAMMPS executable and use Google Test
   v1.6.0. To compile type `make unittests` in the USER-GFMD directory. Then
   execute the "unittests" binary.

2. Compound tests run LAMMPS and check the output. To run these tests
   change to the `tests` subdirectory and run the `run_test.sh` executable with
   the path to the LAMMPS binary as the only argument.
