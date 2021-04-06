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
#ifndef __FORCE_CONSTANTS_H
#define __FORCE_CONSTANTS_H

#include "force.h"
#include "memory.h"

#include "crystal_surface.h"

class ForceConstants {
 public:
    ForceConstants(const char *name, Force *, Error *,
                   const char *pair_style=NULL);

    const char *get_name() const {
        return name_;
    };

    /*!
     * Inform the force constants module about the lattice structure. Allow it
     * to construct supercell.
     */
    virtual void pre_compute(const CrystalSurface *, Error *) { }

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *) = 0;

    /*!
     * On-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
     */
    virtual void onsite(const CrystalSurface *, int, int, double *, Error *);

    /*!
     * Linear force contribution.
     */
    virtual void linear(const CrystalSurface *, double *, Error *error);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *error);

    /*!
     * Return manybody flag
     */
    bool is_manybody() const { return manybody_; };

 protected:
    /*!
     * Descriptor for this surface
     */
    char name_[80];

    /*!
     * Is this a manybody potential? This means that the energy of half of the
     * lattice planes will come from the potential and GFMD, respectively.
     */
    bool manybody_;
};

ForceConstants *force_constants_factory(char *, int, int *, char **,
                                        CrystalSurface *, Force *,
                                        Memory *, Error *);

#endif
