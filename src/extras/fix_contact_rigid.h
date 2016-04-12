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
