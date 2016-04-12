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
#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"

#include "fc_lj_cut_fd.h"

FCLJCutFD::FCLJCutFD(int narg, int *carg, char **arg, CrystalSurface *surface,
                     Force *force, Error *error)
    : FCFiniteDifferences("lj/cut/fd", surface, force, error, "lj/cut/gf")
{
    pair_ = dynamic_cast<LAMMPS_NS::PairLJCutGF*>(force->pair);
    pair_->init();

#ifdef GFMD_DEBUG
    debug_fc_lj_cut_ = new FCLJCut(narg, carg, arg, surface, force, error);
#endif

    surface->compute_neighbor_list(pair_->cut_global);
}


#ifdef GFMD_DEBUG

const double tol = 1e-6;

FCLJCutFD::~FCLJCutFD()
{
    delete debug_fc_lj_cut_;
}


/*!
 * Off-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
 */
void FCLJCutFD::offsite(const CrystalSurface *surf,
                        const CrystalSurface::lattice_neighbor_t *pair,
                        int ab, double *_D, Error *error)
{
    FCFiniteDifferences::offsite(surf, pair, ab, _D, error);

    mat<double> D(3, _D);
    mat<double,3> ref_D = 0.0;
    debug_fc_lj_cut_->offsite(surf, pair, ab, ref_D.data(), error);

    if (pair->ab <= ab) {
        if ((D-ref_D).norm() > tol) {
            std::cout << "offsite, " << * pair << ", ab = " << ab << std::endl;
            D.print();
            ref_D.print();
        }
    }
    else {
        if (D.norm() > tol) {
            std::cout << "offsite, " << * pair << ", ab = " << ab << ", should be zero but is:" << std::endl;
            D.print();
        }
    }
}


/*!
 * On-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
 */
void FCLJCutFD::onsite(const CrystalSurface *surf, int i, int ab,
                       double *_D, Error *error)
{
    FCFiniteDifferences::onsite(surf, i, ab, _D, error);

    mat<double> D(3, _D);
    mat<double,3> ref_D = 0.0;
    debug_fc_lj_cut_->onsite(surf, i, ab, ref_D.data(), error);

    if ((D-ref_D).norm() > tol) {
        std::cout << "onsite, i = " << i << ", ab = " << ab << std::endl;
        D.print();
        ref_D.print();
    }
}


/*!
 * Linear force contribution for the free surface.
 */
void FCLJCutFD::linear(const CrystalSurface *surf, double *linf, Error *error)
{
    FCFiniteDifferences::linear(surf, linf, error);

    vec<double> ref_linf(surf->get_number_of_atoms());
    debug_fc_lj_cut_->linear(surf, ref_linf.data(), error);
    
    if ((ref_linf-linf).nrm2() > tol) {
        vec<double>(surf->get_number_of_atoms(), linf).print();
        ref_linf.print();
    }
}

#endif


/*!
 * Dump generic text info to file
 */
void FCLJCutFD::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    FCFiniteDifferences::dump_info(f, surf, error);
}


/*!
 * Compute forces on all atoms
 */
double FCLJCutFD::energy_and_forces(const CrystalSurface *surf, int *mask,
                                    double *u, double *f, Error *error)
{
    int i,j,jj,jnum,itype,jtype;
    double delx,dely,delz,e,esum,fpair;
    double rsq,r2inv,r6inv,forcelj;

    esum = 0.0;

    const double *x = surf->get_positions();
    const int *type = surf->get_types();
    int nlocal = surf->get_number_of_atoms();

    // loop over neighbors of my atoms

    for (i = 0; i < nlocal; i++) {
        //if (mask && !mask[i]) continue;

        itype = type[i];

        jnum = surf->get_number_of_neighbors(i);
        const CrystalSurface::lattice_neighbor_t *jj = surf->get_neighbors(i);
        for (int njj = 0; njj < jnum; njj++, jj++) {
            j = jj->indexj;
            if (i > j) continue;
            if (mask && !mask[j]) continue;

            delx = jj->r[0]+u[3*i+0]-u[3*j+0];
            dely = jj->r[1]+u[3*i+1]-u[3*j+1];
            delz = jj->r[2]+u[3*i+2]-u[3*j+2];
            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];

            if (rsq < pair_->cutsq[itype][jtype]) {
                r2inv = 1.0/rsq;
                r6inv = r2inv*r2inv*r2inv;
                forcelj = r6inv * (pair_->lj1[itype][jtype]*r6inv - 
                                   pair_->lj2[itype][jtype]);
                fpair = forcelj*r2inv;

                f[3*i+0] += delx*fpair;
                f[3*i+1] += dely*fpair;
                f[3*i+2] += delz*fpair;
                f[3*j+0] -= delx*fpair;
                f[3*j+1] -= dely*fpair;
                f[3*j+2] -= delz*fpair;

                e = r6inv*(pair_->lj3[itype][jtype]*r6inv-
                           pair_->lj4[itype][jtype]) - 
                    pair_->offset[itype][jtype];
                if (i == j)  e /= 2;
                esum += e;
            }
        }
    }

    return esum;
}
