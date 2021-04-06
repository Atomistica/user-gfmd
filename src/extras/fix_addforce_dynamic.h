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

FixStyle(addforce/dynamic,FixAddForceDynamic)

#else

#ifndef LMP_FIX_ADDFORCE_DYNAMIC_H
#define LMP_FIX_ADDFORCE_DYNAMIC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddForceDynamic : public Fix {
 public:
  FixAddForceDynamic(class LAMMPS *, int, char **);
  ~FixAddForceDynamic();
  int setmask();
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double compute_vector(int);

 private:
  double xvalue_, yvalue_, zvalue_;   // force amplitude and direction
  int nx_, ny_;                       // Fourier components
  double freq_;                       // driving frequency

  int nsteps_;                        // number of steps

  int vec_is_reduced_;
  double vec_[12], reduced_vec_[12];  // Fourier components of position and vel.
};

}

#endif
#endif
