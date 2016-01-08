#ifndef LMP_COMM_H
#define LMP_COMM_H

namespace LAMMPS_NS {

class Comm {
 public:
  int me,nprocs;                    // proc info
  int procgrid[3];                  // assigned # of procs in each dim

  Comm()
  {
    me = 0;
    nprocs = 1;
    procgrid[0] = procgrid[1] = procgrid[2] = 1;
  }
};

}

#endif
