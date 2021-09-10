/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(gfmd,FixGFMD)

#else

#ifndef FIX_GFMD_H
#define FIX_GFMD_H

#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "fix.h"
#include "linearalgebra.h"

#include "gfmd_solver.h"
#include "surface_stiffness.h"
	 
namespace LAMMPS_NS {

class FixGFMD : public Fix {
 
 public:
  FixGFMD(class LAMMPS *, int, char **);
  virtual ~FixGFMD();

  virtual int  setmask();
  virtual void init();
  virtual void setup(int);
  virtual void min_setup(int);
  virtual void pre_force(int);
  virtual void post_force(int);
  virtual void min_pre_force(int);
  virtual void min_post_force(int);
  virtual double memory_usage();
  virtual double compute_scalar();
  virtual double compute_vector(int);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *); 

 protected:
  int me, nprocs;                // my process rank, total number of proc (MPI)
  int origin_tag;                // atomic tag of the surface origin for FFT

  bool do_dump_stiffness;        // dump the stiffness coefficients to files?
  bool do_dump_greens_function;  // dump the greens function to files?
  int dumpr_every;               // dump interval for displ., for. energies
  int dumpq_every;               // dump interval for displ., for. energies

  double masstotal, rmasstotal;  // total mass of GF atoms and its inverse

  int noutfor;                   // frequency to output Green's Function force on GF atoms;
  char *prefix;                  // prefix for output files in GFMD
  FILE *gfmdlog;                 // file to output GF related info

  double fsum_loc[3];            // to store the original forces
  double fsum[3];
  bool   fsum_flag;              // indicates if all forces have been collected or not
  double epot;                   // effective potential energy of GFMD layer
  bool   pre_flag;               // has pre_force been called?
  bool   epot_flag;              // has epot been computed?

  bool   grid_checked_flag;      // have grid indices been checked?

  int ndof, ndof_sq;             // kernel dimension

  int nx, ny, nu, rnu;           // dimension of FFT mesh
  int nxy;                       // nx*ny
  int natoms;                    // total number of GF atoms

  int nmax;                      // maximum number of local atoms

  int xlo_loc, xhi_loc;          // part of grid on local proc, x-dim
  int ylo_loc, yhi_loc;          // part of grid on local proc, y-dim
  int nx_loc, ny_loc;            // number of grid points in x-, y-dim
  int nxy_loc;                   // total number of grid points on this proc
  int natoms_localsubgrid;       // number of GF atoms locally

  int ngfmd_loc;                 // counter for GFMD atoms on curr. domain

  int xshift, yshift;            // shift of center of mass in lattice coord. 
  //int low_left_layer;          // GF layer, counting from top, which contains the
                                 // atom in each unit cell that is most "lower left-hand"
                                 // (user inputs this to match their geometry file.  default is 
                                 //  nu-1, typically the minimum-z layer.) *not implemented as not needed now*

  int mynq, mynpt;

  bool fix_center[3];            // fix center of mass?

  bool reset_xeq;                // flag to indicate whether/how to reset xeq according to surface lattice info from gfc
  bool reset_map;

  double xeqcm[3];               // center of mass of xeq

  double pre_time, post_time;    // timings
  bigint last_timestep;          // timestep of last GFMD call

  double border;                 // border for communication

  double **u_xy;                 // displacement field, as a function of xy
  double **f_xy;                 // force field, as a function of xy index
  double **f_i;                  // force field, as a function of atom index

  StiffnessKernel *kernel;       // stiffness kernel
  GFMDSolver *solver;            // solver instance

  bool *force_is_known;          // did we compute or receive this force?

  void grid_to_list();

  // method to compute the initial configuration
  void comp_xeq();
  // method to compute mapping info from initial configuration based on
  // surface vectors
  void comp_map();
  void check_grid_indices();

  void get_cm_from_xeq(double *);

  void dump();
};

}

#endif
#endif
