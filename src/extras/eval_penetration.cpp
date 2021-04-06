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
/*
 * Standalone code, count number of atoms that penetrate a contact map.
 */
#include <stdio.h>

#include <netcdfcpp.h>

#include "error.h"
#include "memory.h"
#include "table2d.h"

using namespace std;
using namespace LAMMPS_NS;

void read_map(char *fn, double z0, int &nx, int &ny, double **&map_data,
	      Memory *memory, Error *error)
{
  FILE *f = fopen(fn, "r");

  if (!f) {
    char errstr[1024];
    sprintf(errstr, "Could not open contact map file '%s'.", fn);
    error->all(errstr);
  }

  // count number of columns
  nx = 0;
  char c = fgetc(f), lastc = ' ';
  while (!feof(f) && c != '\n') {
    if ((isalnum(c) || c == '-' || c == 'e') && isblank(lastc))
      nx++;
    lastc = c;
    c = fgetc(f);
  }
  if (feof(f))
    error->all("fix contact/map: Could not read map file.");

  if (nx <= 0)
    error->all("fix contact/map: Something wrong: nx <= 0");


  /*
   * for now assume a square map
   */
  ny = nx;

  rewind(f);

  map_data = memory->create_2d_double_array(nx, ny, "map");
  double za = 0.0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      double tmp;
      if (fscanf(f, "%lf", &tmp) != 1)
	error->all("fix contact/map: Error read map file.");

      map_data[i][j] = tmp;

      za += tmp;
    }
  }

  fclose(f);

  // remove average height and add zero height
  za = z0 - za/(nx*ny);
  for (int i = 0; i < nx*ny; i++) {
    map_data[0][i] += za;
  }
}


int main(int narg, char **argv)
{
  Error *error = new Error();
  Memory *memory = new Memory(error);

  double z0;
  char *endptr;
  
  z0 = strtod(argv[3], &endptr);
  if (endptr == argv[3])
    error->all("Could not convert z0 argument to double.");

  int nx, ny;
  double **map_data;
  read_map(argv[2], z0, nx, ny, map_data, memory, error);
  Table2D *map = new Table2D(nx, ny, map_data, true, error, memory);
  memory->destroy_2d_double_array(map_data);

  NcFile *nc = new NcFile(argv[1]);

  int nframes = nc->get_dim("frame")->size();
  int natoms = nc->get_dim("atom")->size();

  NcVar *x_var = nc->get_var("coordinates");
  double **x = memory->create_2d_double_array(natoms, 3, "x");
  double *x_ptr = x[0];

  NcVar *l_var = nc->get_var("cell_lengths");


  for (int frame = 0; frame < nframes; frame++) {
    if (!x_var->set_cur(frame))
      error->all("set_cur failed.");
    if (!x_var->get(x_ptr, 1, natoms, 3))
      error->all("get failed.");

    double l[3];
    if (!l_var->set_cur(frame))
      error->all("set_cur failed.");
    if (!l_var->get(l, 1, 3))
      error->all("get failed.");

    int n = 0;
    for (int i = 0; i < natoms; i++) {
      double rx, ry;
      double z, dz_dx, dz_dy;

      rx = x[i][0]*nx/l[0];
      ry = x[i][1]*ny/l[1];

      map->eval(rx, ry, z, dz_dx, dz_dy);

      if (x[i][2] > z)
	n++;
    }

    printf("%i %i %f\n", frame, n, double(n)/natoms);
  }

  memory->destroy_2d_double_array(x);

  delete nc;
  //  delete map;
  delete memory;
  delete error;
}
