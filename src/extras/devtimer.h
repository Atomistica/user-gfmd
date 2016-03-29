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

#ifndef LMP_DEVTIMER_H
#define LMP_DEVTIMER_H

#include "pointers.h"

// Add timing fields as desired to the start of this enumerated list.
enum{DEVTIME_GFINIT,
     DEVTIME_GFSETUP,
     DEVTIME_GFPREFORCE,
     DEVTIME_GFPOSTFORCE,
     DEVTIME_TOTAL,
     DEVTIME_GETD,
     DEVTIME_GETG,
     DEVTIME_TIMING,
     DEVTIME_PREFORCE,
     DEVTIME_POSTFORCE,
     DEVTIME_N};

namespace LAMMPS_NS {

class Devtimer : protected Pointers {
 public:
  double *devtimetotals;
  double *devtimestarts;

  Devtimer(class LAMMPS *);
  ~Devtimer();
  void start(int);
  void stop(int);
  void reset(int);

};

}

#endif
