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
#ifndef __FT_STIFFNESS_H
#define __FT_STIFFNESS_H

#include "crystal_surface.h"
#include "domain.h"
#include "error.h"
#include "force_constants.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class FTStiffnessKernel : public StiffnessKernel {
 public:
    FTStiffnessKernel(int, int *, char **, Domain *, Force *, Memory *,
                      Error *);
    ~FTStiffnessKernel();

    virtual const double *get_cell() const {
        return crystal_surface_->get_cell();
    }

    virtual const double *get_invcell() const {
        return crystal_surface_->get_invcell();
    }

    virtual bool get_mapping(int num_ref, double *ref_positions, int *indices,
                             int *types=NULL, double tol=1e-3) const {
        return crystal_surface_->get_mapping(num_ref, ref_positions, indices,
                                             types, tol);
    }

    /*!
     * Called before get_dynamical_matrices and used to tabulate Dij for each
     * pair ij.
     */
    virtual void pre_compute();

    /*!
     * Get per cell force-constant matrices U and V.
     */
    virtual void get_dynamical_matrices(double, double, double_complex *,
                                        double_complex *, double_complex *,
                                        double_complex dU=0.0);

    virtual void get_force_at_gamma_point(double *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *);

    /*!
     * Is this a manybody potential? This means that the energy of half of the
     * lattice planes will come from the potential and GFMD, respectively.
     */
    virtual bool is_manybody() { return force_constants_->is_manybody(); }

 private:
    /*
     * Surface normal
     */
    int z_;

    /*
     * Forces added to each surface layer
     * The per-atom contribution, normal to the surface, to
     * each layer to offset surface relaxation.  Max interaction range
     * implemented is 2 atomic layers, so max many-body nu_ is 2*2.
     */
    double linf_[2*2];

    /*!
     * Crystal lattice
     */
    CrystalSurface *crystal_surface_;

    /*!
     * Interation model, force-constant matrix
     */
    ForceConstants *force_constants_;

    /*!
     * The actual values of the force constants between pairs in the interfacial
     * layer, on-site.
     */
    double *D0ii_;

    /*!
     * The actual values of the force constants between pairs in the bulk,
     * on-site.
     */
    double *D1ii_;

    /*!
     * The actual values of the force constants between pairs in the interfacial
     * layer, off-site.
     */
    double *D0ij_;

    /*!
     * The actual values of the force constants between pairs in the bulk,
     * off-site.
     */
    double *D1ij_;
};

}

#endif
