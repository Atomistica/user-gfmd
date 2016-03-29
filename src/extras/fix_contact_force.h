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

FixStyle(contact/map,FixContactForce)

#else

#ifndef LMP_FIX_CONTACT_FORCE_H
#define LMP_FIX_CONTACT_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactForce : public Fix {
 public:
  FixContactForce(class LAMMPS *, int, char **);
  ~FixContactForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void min_pre_force(int);

  double compute_scalar();
  double compute_vector(int);

 private:
  double **forces;                // external forces
 
  int force_flag;                 // have the forces been comm.?
  double force[4], force_loc[4];  // total force on rigid body

  void read_forces(char *);
};

}

#endif
#endif
