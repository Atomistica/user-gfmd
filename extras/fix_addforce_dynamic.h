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
