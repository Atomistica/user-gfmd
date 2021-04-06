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
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contact/sphere,FixContactSphere)

#else

#ifndef LMP_FIX_CONTACT_SPHERE_H
#define LMP_FIX_CONTACT_SPHERE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactSphere : public Fix {
 public:
  FixContactSphere(class LAMMPS *, int, char **);
  ~FixContactSphere();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);

  double compute_scalar();
  double compute_vector(int);

 private:
  char *cxvar, *cyvar, *czvar;
  double cx, cy, cz, radius, epsilon, sigma, cutoff;  // parameters
  double radius_sq, radius_plus_cutoff_sq;
  double coeff1, coeff2, coeff3, coeff4, offset;

  int force_flag;                            // have the forces been comm.?
  double fsphere[5], fsphere_loc[5];         // total force of sphere body
};

}

#endif
#endif
