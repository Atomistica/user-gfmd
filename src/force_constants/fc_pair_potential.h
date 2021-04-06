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

FCStyle(pair-potential,FCPairPotential)

#else

#ifndef __FC_PAIR_POTENTIAL_H
#define __FC_PAIR_POTENTIAL_H

#include "force_constants.h"

class FCPairPotential : public ForceConstants {
 public:
    FCPairPotential(int, int *, char **, CrystalSurface *, Force *, Error *);
    virtual ~FCPairPotential();

    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *);

 protected:
    /*!
     * Spring constants for nearest- and next-nearest neighbors
     * k is d2V/dr^2_ij at equilibrium r_ij  (Eq A4)
     * kappa is dV/dr_ij * 1/r_ij at equilibrium r_ij (Eq A6)
     */
    double k_[3], kappa_[3];

    /*!
     * Linear force contribution (for each layer)
     */
    double *linf_;
};

#endif

#endif
