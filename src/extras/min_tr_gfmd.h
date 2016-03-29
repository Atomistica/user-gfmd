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

#ifdef MINIMIZE_CLASS

MinimizeStyle(tr/gfmd,MinTRGFMD)

#else

#ifndef LMP_MIN_TR_CONTACT_H
#define LMP_MIN_TR_CONTACT_H

#include "min_linesearch.h"

namespace LAMMPS_NS {

class MinTRGFMD : public MinLineSearch {
 public:
  MinTRGFMD(class LAMMPS *);
  void setup_style();
  int iterate(int);

 private:
  class FixGFMD *gfmd;
};

}

#endif
#endif
