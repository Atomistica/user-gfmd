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
#ifndef LMP_NC_TRAJ_IO_H
#define LMP_NC_TRAJ_IO_H

#include "pointers.h"

namespace LAMMPS_NS {

const int NC_FIELD_NAME_MAX = 100;

class NCTrajIO : protected Pointers {
 public:
  NCTrajIO(class Pointers *);
  ~NCTrajIO();

  void open_for_read(char *);
  void read_frame(int);

  void open_for_write(char *);
  void write_header(double, double, double, double);
  void write_data();

 private:
  int framei;                  // current frame index

  int ntotal;                  // # of atoms

  bool double_precision;       // write everything as double precision

  int n_buffer;                // size of buffer
  int *int_buffer;             // buffer for passing data to netcdf
  double *double_buffer;       // buffer for passing data to netcdf

  int ncid;

  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;
  int cell_angular_dim;
  int label_dim;

  int spatial_var;
  int cell_spatial_var;
  int cell_angular_var;

  int time_var;
  int cell_lengths_var;
  int cell_angles_var;
  int epot_var;
  int evdwl_var;
  int eel_var;
  int eel0_var;

  int Z_var;
  int coordinates_var;
  int forces_var;

  void ncerr(int, int);
};

}

#endif
