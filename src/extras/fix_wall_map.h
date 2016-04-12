/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

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

#ifdef FIX_CLASS

FixStyle(wall/map,FixWallMap)

#else

#ifndef LMP_FIX_WALL_MAP_H
#define LMP_FIX_WALL_MAP_H

#include "fix.h"
#include "table2d.h"

namespace LAMMPS_NS {

class FixWallMap : public Fix {
 public:
  FixWallMap(class LAMMPS *, int, char **);
  ~FixWallMap();
  int setmask();
#ifndef NO_LAMMPS
  void init();
  void setup(int);
#endif
  void min_setup(int);
  void pre_force(int);
  void min_pre_force(int);

  void write_restart(FILE *);
  void restart(char *);

  double compute_scalar();
  double compute_vector(int);

  double get_rms_height() {
    return h0_;
  }

  double curvt[6];               // curvature: xx, yy, zz, xy, yz, xz

 private:
  // potential parameters
  double epsilon_, sigma_, cutoff_, inner_cutoff_, mcutoff_; 

  double mc0_, mc1_, mc2_;
  double coeff1_, coeff2_, coeff3_, coeff4_, coeff5_, coeff6_, offset_;
  double smcff0_, smcff1_, smcff2_, smcff3_, smcff4_;

  double lambda_;                // short wavelength cutoff

  bool firstcall_;               // is this the first call?
  double pendist_;               // reset atom if surface is penetrated

  double z0_;                    // average height
  int z0_type_;                  // is z0 average or min position

  double x0_, y0_, vx0_, vy0_;   // map offset and its velocity

  int norm, nx_, ny_;            // surface normal, map dimensions
  int xlo_loc0_, xhi_loc0_;      // part of grid on local proc, x-dim
  int ylo_loc0_, yhi_loc0_;      // part of grid on local proc, y-dim
  int nx_loc0_, ny_loc0_;        // number of grid points in x-, y-dim
  int xlo_loc_, xhi_loc_;        // part of grid on local proc, x-dim, with buf
  int ylo_loc_, yhi_loc_;        // part of grid on local proc, y-dim, with buf
  int nx_loc_, ny_loc_;          // number of grid points in x-, y-dim, buf
  int nxy_loc0_, nxy_loc_;       // total number of grid points on this proc

  double **map_data_;            // map data
  Table2D *map_;                 // map object

  double h0_;                    // rms height for maps, amplitude for cryst.
  double qx_, qy_;               // reciprocal lattice vectors

  int force_flag;                // have the forces been comm.?
  // this is: e, fx, fy, fz, npen, nrep, natt
  double fmap_[7], fmap_loc_[7];

  void (FixWallMap::*eval_pot_func_)(double, double &, double &, double &);

  void init_epsilon_and_sigma(int &, int, char **);
  void init_cutoff(int &, int, char **);
  void init_outer_cutoff(int &, int, char **);

  void init_lj93(int &, int, char **);
  void eval_lj93(double, double &, double &, double &);

  void init_lj126(int &, int, char **);
  void eval_lj126(double, double &, double &, double &);

  void init_lj93smooth(int &, int, char **);
  void eval_lj93smooth(double, double &, double &, double &);

  void init_exp(int &, int, char **);
  void eval_exp(double, double &, double &, double &);

  void init_expf(int &, int, char **);
  void eval_expf(double, double &, double &, double &);

  void init_dugdale(int &, int, char **);
  void eval_dugdale(double, double &, double &, double &);

  void (FixWallMap::*eval_map_func_)(double, double, double &, double &,
                                     double &, double &, double &, double &);

  void init_flat();
  void eval_flat(double, double, double &, double &, double &, double &,
                 double &, double &);

  void init_100(int &, int, char **);
  void eval_100(double, double, double &, double &, double &, double &,
                double &, double &);

  void init_111(int &, int, char **);
  void eval_111(double, double, double &, double &, double &, double &,
                double &, double &);

  void allocate_map(double **&, int, int);
  void read_map(char *, double **&);
  void write_map(char *, double **);
  void filter_map(double **);
  void displace_map(double **);
  void eval_table2d(double, double, double &, double &, double &, double &,
                    double &, double &);
  void eval_raw_data(double, double, double &, double &, double &, double &,
                     double &, double &);
};

}

#endif
#endif
