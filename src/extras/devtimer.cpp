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

// A basic stopwatch that developers can use to time parts of the code.
// Start, stop, and reset can be inserted into the c++ source code.
// All that is required is the devtimer.h include at the file head.
// As long as class constructor has : Pointers(lmp) and lmp is an input.
// Can be printed out with the normal timing in finish.cpp.
// In the long term, this functionality could be built into LAMMPS timer.
// For now, a separate devtimer is easier to keep integrating into LAMMPS
// Tristan Sharp 7/26/12

#include "mpi.h"
#include "devtimer.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Devtimer::Devtimer(LAMMPS *lmp) : Pointers(lmp)
{
  memory->create(devtimetotals,DEVTIME_N,"array");
  memory->create(devtimestarts,DEVTIME_N,"array");
  for (int i = 0; i < DEVTIME_N; i++) devtimetotals[i] = 0.0;  
  for (int i = 0; i < DEVTIME_N; i++) devtimestarts[i] = 0.0;  
}



Devtimer::~Devtimer()
{
  memory->destroy(devtimetotals);
  memory->destroy(devtimestarts);
}
// MPI_Barrier(world) could sync procs timing, but why?
void Devtimer::start(int which)
{  devtimestarts[which] = MPI_Wtime();}

void Devtimer::stop(int which)
{  devtimetotals[which] += (double)(MPI_Wtime()) - devtimestarts[which];}

void Devtimer::reset(int which)
{  devtimetotals[which] = 0.0;}



