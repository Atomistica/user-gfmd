#ifndef DOMAIN_H
#define DOMAIN_H

namespace LAMMPS_NS {

class Domain {
 public:
  int triclinic;		         // 0 = orthog box, 1 = triclinic

  double xprd,yprd,zprd;                 // global box dimensions
  double boxlo[3],boxhi[3];              // orthogonal box global bounds
  double sublo[3],subhi[3];              // sub-box bounds on this proc

  Domain()
  {
    triclinic = false;
    set_cell(1,1,1);
  }

  void set_cell(double x, double y, double z)
  {
    xprd = x;
    yprd = y;
    zprd = z;

    boxlo[0] = boxlo[1] = boxlo[2] = 0;
    boxhi[0] = x;
    boxhi[1] = y;
    boxhi[2] = z;

    sublo[0] = sublo[1] = sublo[2] = 0;
    subhi[0] = x;
    subhi[1] = y;
    subhi[2] = z;
  }
};

};

#endif
