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

#ifndef LMP_STUB_ERROR_H
#define LMP_STUB_ERROR_H

#include <stdio.h>
#include <stdlib.h>

#define FLERR __FILE__,__LINE__

namespace LAMMPS_NS {

class Error {
 public:
  void all(const char *file, int line, const char *errstr) {
    printf("all: %s\n", errstr);
    exit(999);
  }

  void one(const char *file, int line, const char *errstr) {
    printf("one: %s\n", errstr);
    exit(999);
  }
};

}

#endif
