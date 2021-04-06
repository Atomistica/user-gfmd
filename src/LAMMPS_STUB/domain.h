/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

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
