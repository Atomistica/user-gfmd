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
