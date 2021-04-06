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
#ifdef FC_CLASS

FCStyle(lj/cut,FCLJCut)

#else

#ifndef __FC_LJ_CUT_H
#define __FC_LJ_CUT_H

#include "force_constants.h"

#include "pair_lj_cut_gf.h"

class FCLJCut : public ForceConstants {
 public:
    FCLJCut(int, int *, char **, CrystalSurface *, Force *, Error *);

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *, Error *error);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * First derivative of pair expression
     */
    virtual double first_derivative(int, int, double);

    /*!
     * Second derivative of pair expression
     */
    virtual double second_derivative(int, int, double);

    /*!
     * The lj/cut pair style class
     */
    LAMMPS_NS::PairLJCutGF *pair_;
};

#endif

#endif
