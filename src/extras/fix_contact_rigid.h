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

FixStyle(contact/rigid,FixContactRigid)

#else

#ifndef LMP_FIX_CONTACT_RIGID_H
#define LMP_FIX_CONTACT_RIGID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactRigid : public Fix {
 public:
  FixContactRigid(class LAMMPS *, int, char **);
  ~FixContactRigid();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);

  double compute_scalar();
  double compute_vector(int);

  double memory_usage();

 private:
  int output_interval, nit, nout;
  bool output_forces;                      // write forces to file
  char *prefix;                            // file prefix

  int rank;                                // rank of this processor

  int np;                                  // number of p. for pressure hist.
  double minp, maxp;                       // minimum, maximum pressure
  double minpsq, maxpsq;

  int ncontact_loc, ncontact_sum;          // number of atoms in contact
  int force_flag;                          // have the forces been comm.?
  double frigid[3], frigid_all[3];         // total force of rigid body

  double *phist_loc, *phist_sum;           // pressure histogram
};

}

#endif
#endif
