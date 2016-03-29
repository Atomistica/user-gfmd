#ifndef LMP_ATOM_H
#define LMP_ATOM_H

#include "lmptype.h"
#include "memory.h"

namespace LAMMPS_NS {

class Atom {
 public:
  bigint natoms;                // total # of atoms in system, could be 0
  int nlocal,nghost;            // # of owned and ghost atoms on this proc
  int nmax;                     // max # of owned+ghost in arrays on this proc

  // per-atom arrays
  // customize by adding new array

  int *tag,*type,*mask;
  double **x,**v,**f,**xeq;

  Atom(bigint, Memory *);
  virtual ~Atom();

 protected:
  Memory *memory;
};

}

#endif
