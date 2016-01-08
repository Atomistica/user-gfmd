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

FixStyle(gfmd/static,FixGFMDStatic)

#else

#ifndef FIX_GFMD_STATIC_H
#define FIX_GFMD_STATIC_H

#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "fix.h"
#include "linearalgebra.h"

#include "gfmd_solver.h"
#include "gfmd_stiffness.h"
	 
namespace LAMMPS_NS {

class FixGFMDStatic : public Fix {
 
 public:
  FixGFMDStatic(class LAMMPS *, int, char **);
  ~FixGFMDStatic();

  int  setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void min_pre_force(int);
  void min_post_force(int);
  void end_of_step();
  double memory_usage();
  double compute_scalar();
  double compute_vector(int);

 private:
  int me, nprocs;                // my process rank, total number of proc (MPI)
  int origin_tag;                // atomic tag of the surface origin for FFT

  bool do_dump_stiffness;        // dump the stiffness coefficients to files?

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

  int sysdim, surfdim;           // system dimension (e.g. 3) and
                                 // surface dimension (e.g. 2)
  int ndof, ndof_sq;             // kernel dimension

  int nx, ny, nu;                // dimension of FFT mesh
  int nxy;                       // nx*ny
  int natoms;                    // total number of GF atoms

  int nmax;                      // maximum number of local atoms

  int xlo_loc, xhi_loc;          // part of grid on local proc, x-dim
  int ylo_loc, yhi_loc;          // part of grid on local proc, y-dim
  int nx_loc, ny_loc;            // number of grid points in x-, y-dim
  int nxy_loc;                   // total number of grid points on this proc
  int natoms_loc;                // number of GF atoms locally
  int natoms_loc_recv;           // number of local GF atoms received

  int ngfmd_loc;                 // counter for GFMD atoms on curr. domain

  int xshift, yshift;            // shift of center of mass in lattice coord.

  int mynq, mynpt;

  bool have_su, have_sv;         // have the surface vectors been specified?
  double surfvec[2][2];          // surface vectors

  bool fix_center[3];            // fix center of mass?

  bool reset_xeq;                // flag to indicate whether/how to reset xeq according to surface lattice info from gfc
  bool reset_map;

  double xeqcm[3];               // center of mass of xeq

  double pre_time, post_time;    // timings

  double **u_xy;                 // displacement field
  double **f_xy;                 // force field

  double_complex **phi_q;        // stiffness kernel buffer

  GFMDSolver *solver;            // solver instance

  bool *force_is_known;          // did we compute are receive this force?

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *); 

  // Phi related methods and variables
  void readphi(char *);

  // method to compute the initial configuration
  void comp_xeq();
  // method to compute mapping info from initial configuration based on
  // surface vectors
  void comp_map(int, int *, int *);

  void phi_q0_rescale(double_complex *);
  void phi_analytic(StiffnessKernel *kernel);

  void dump_stiffness();

  void get_cm_from_xeq(double *);
};

}

#endif
#endif
