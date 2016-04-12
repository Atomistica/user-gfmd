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
#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"

#include "fc_lj_cut.h"

FCLJCut::FCLJCut(int narg, int *carg, char **arg, CrystalSurface *surface,
                 Force *force, Error *error)
    : ForceConstants("lj/cut", force, error, "lj/cut/gf")
{
    pair_ = dynamic_cast<LAMMPS_NS::PairLJCutGF*>(force->pair);
    pair_->init();

    surface->compute_neighbor_list(pair_->cut_global);
}


/*!
 * Dump generic text info to file
 */
void FCLJCut::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    ForceConstants::dump_info(f, surf, error);

    int nat = surf->get_number_of_atoms();

    fprintf(f, "Force constants for each neighbor shell:\n");
    for (int i = 0; i < nat; i++) {
        int itype = surf->get_type(i);
        fprintf(f, "Atom %i, k (kappa) = ", i);
        for (int j = 0; j < 10; j++) {
            const CrystalSurface::lattice_neighbor_t *neigh = 
                surf->get_neighbors(i);
            int nneigh = surf->get_number_of_neighbors(i);
            int n = 0;
            while (n < nneigh && neigh->neighbor_shell != j) {
                n++;
                neigh++;
            }
            if (n < nneigh) {
                int jtype = surf->get_type(neigh->indexj);
                fprintf(f, "%f (%f) ",
                        second_derivative(itype, jtype, neigh->rnorm),
                        first_derivative(itype, jtype, neigh->rnorm)/
                        neigh->rnorm);
            }
        }
        fprintf(f, "\n");
    }

    double linf[nat];
    linear(surf, linf, error);

    fprintf(f, "Linear forces on each atom = ");
    for (int i = 0; i < nat; i++)  fprintf(f, "%f ", linf[i]);
    fprintf(f, "\n");
}


/*!
 * Off-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
 */
void FCLJCut::offsite(const CrystalSurface *surf,
                      const CrystalSurface::lattice_neighbor_t *pair,
                      int ab, double *_D, Error *error)
{
    int itype = surf->get_type(pair->indexi);
    int jtype = surf->get_type(pair->indexj);

    vec<double, 3> e(pair->r);
    e /= e.nrm2();

    mat<double> D(3, _D);
    outer(e, e, D);

    mat<double> perp = identity<double>(3) - D;

    double k = second_derivative(itype, jtype, pair->rnorm);
    double kappa = first_derivative(itype, jtype, pair->rnorm)/pair->rnorm;

    D *= -k;
    perp *= kappa;

    D -= perp;
}


/*!
 * First derivative of pair expression
 */
double FCLJCut::first_derivative(int itype, int jtype, double r)
{
    double rsq = r*r;
    double r2inv = 1.0/rsq;
    if (rsq < pair_->cut[itype][jtype]*pair_->cut[itype][jtype]) {
        double r6inv = r2inv*r2inv*r2inv;
        return r6inv * (pair_->lj2[itype][jtype] -
                        pair_->lj1[itype][jtype]*r6inv) / r;
    }
    return 0.0;
}


/*!
 * Second derivative of pair expression
 */
double FCLJCut::second_derivative(int itype, int jtype, double r)
{
    double rsq = r*r;
    double r2inv = 1.0/rsq;
    if (rsq < pair_->cut[itype][jtype]*pair_->cut[itype][jtype]) {
        double r6inv = r2inv*r2inv*r2inv;
        return r6inv * r2inv * (13*pair_->lj1[itype][jtype]*r6inv-
                                7*pair_->lj2[itype][jtype]);
    }
    return 0.0;
}


/*!
 * Linear force contribution for the free surface.
 */
void FCLJCut::linear(const CrystalSurface *surf, double *linf, Error *error)
{
    mat<double,3> cell(surf->get_cell());

    int nat = surf->get_number_of_atoms();
    std::fill(linf, linf+nat, 0.0);

    for (int i = 0; i < nat; i++) {
        double f = 0.0;
        int itype = surf->get_type(i);
        const CrystalSurface::lattice_neighbor_t *neigh = 
            surf->get_neighbors(i);
        for (int n = 0; n < surf->get_number_of_neighbors(i); n++) {
            if (neigh->ab <= 0) {
                int jtype = surf->get_type(neigh->indexj);
                vec<double,3> e = neigh->r;
                e /= neigh->rnorm;
                f += first_derivative(itype, jtype, neigh->rnorm)*e[2];
            }
            neigh++;
        }
        linf[i] = -f;
    }
}
