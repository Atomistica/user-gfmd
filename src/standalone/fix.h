#ifndef FIX_H
#define FIX_H

#include "pointers.h"

namespace LAMMPS_NS {

class Fix : protected Pointers {
 public:
  int groupbit;

  int nevery;                    // how often to call an end_of_step fix

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int size_vector;               // length of global vector
  int global_freq;               // frequency s/v data is available at

  int extscalar;            // 0/1 if global scalar is intensive/extensive
  int extvector;            // 0/1/-1 if global vector is all int/ext/extlist

  int INITIAL_INTEGRATE,POST_INTEGRATE;    // mask settings
  int PRE_EXCHANGE,PRE_NEIGHBOR;
  int PRE_FORCE,POST_FORCE,FINAL_INTEGRATE,END_OF_STEP,THERMO_ENERGY;
  int INITIAL_INTEGRATE_RESPA,POST_INTEGRATE_RESPA;
  int PRE_FORCE_RESPA,POST_FORCE_RESPA,FINAL_INTEGRATE_RESPA;
  int MIN_PRE_EXCHANGE,MIN_PRE_FORCE,MIN_POST_FORCE,MIN_ENERGY;
  int POST_RUN;

  Fix(class Pointers *ptrs, int, char **): Pointers(ptrs) {
    INITIAL_INTEGRATE = 1;
    POST_INTEGRATE = 2;
    PRE_EXCHANGE = 4;
    PRE_NEIGHBOR = 8;
    PRE_FORCE = 16;
    POST_FORCE = 32;
    FINAL_INTEGRATE = 64;
    END_OF_STEP = 128;
    THERMO_ENERGY = 256;
    INITIAL_INTEGRATE_RESPA = 512;
    POST_INTEGRATE_RESPA = 1024;
    PRE_FORCE_RESPA = 2048;
    POST_FORCE_RESPA = 4096;
    FINAL_INTEGRATE_RESPA = 8192;
    MIN_PRE_EXCHANGE = 16384;
    MIN_PRE_FORCE = 32768;
    MIN_POST_FORCE = 65536;
    MIN_ENERGY = 131072;
    POST_RUN = 262144;

    groupbit = 1;
  }
  virtual ~Fix() { }
};

}

#endif
