Green's function molecular dynamics code
========================================

If you use this code, please cite:

- C. Campana, M.H. Muser, Phys. Rev. B 74, 075420 (2006)
- L. Pastewka, T.A. Sharp, M.O. Robbins, Phys. Rev. B 86, 075459 (2012)

LAMMPS
------

The files in this directory are a user-contributed package for LAMMPS.

This package implements the Green's function molecular dynamics, by
introducing one more fixes to LAMMPS. `fix gfmd`, is used to evaluate the
 effective elastic forces by using the elastic stiffness coefficients.

See the documentation files of this commands for detailed usage.

The compiling of this package along with LAMMPS requires the FFT3d wrappers
from the kspace package of LAMMPS be included as well. 
An example of the Makefile is also available in this directory.

EAM GF capability requires "make yes-MANYBODY", which places 
pair_eam.cpp and other files into the src/ directory, and then 
make yes-USER-GFMD will patch pair_eam.cpp.

Besides, you should also append 'user-gfmd' to the 'PACKUSER' variable
in the file of src/Makefile before invoke `make yes-user`. Note that
`make yes-user-gmfd` does not work if used with the USER-CUDA package.
If used with USER-CUDA, go into the USER-GFMD directory and execute
`sh Install.sh 1` manually.

Functionality beyond pair potentials will be included soon.

Testing
-------

Parts of the code are included in a testing framework that uses Google Test
v1.6.0. To compile type

make unittests

in the USER-GFMD directory and execute the "unittests" binary.

Version compatibility
---------------------

- Up to rev 161 works with lammps-5Jan11
- Up to rev XXX works with lammps-14Aug11
