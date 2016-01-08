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
  void setup(int); // TAS added 4/17/14 to see if can dump fix forces
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 protected:
  double shiftx_nn_; // shift of hypothetical opposing lattice, units nn
  double shifty_nn_; // shift of hypothetical opposing lattice, units nn
  double acontact_nn_;  // effective contact radius, units of nn spacing
  double nndist_;
  double r_overEstar;

  int force_flag;
  double foriginal_all[4];
  double foriginal[4];
};

}

#endif
#endif
