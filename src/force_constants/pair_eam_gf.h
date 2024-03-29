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

#ifdef PAIR_CLASS

PairStyle(eam/gf,PairEAMGF)

#else

#ifndef LMP_PAIR_EAM_GF_H
#define LMP_PAIR_EAM_GF_H

#include <stdio.h>
#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMGF : public PairEAM {
 public:

  PairEAMGF(class LAMMPS *);
  virtual ~PairEAMGF();
  virtual void compute(int, int);
};

}

#endif
#endif

