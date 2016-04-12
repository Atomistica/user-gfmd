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
// Eqn labels refer to Pastewka et. al. PRB 86 075459 (2012)

#define _USE_MATH_DEFINES
#include <iostream>

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "ft_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Face centered cubic (100) surface
 * --------------------------------------------------------------------*/

FTStiffnessKernel::FTStiffnessKernel(int narg, int *carg, char **arg,
                                     Domain *domain, Force *force,
                                     Memory *memory, Error *error)
    : StiffnessKernel(narg, carg, arg, domain, memory, error), D0ii_(NULL),
      D1ii_(NULL), D0ij_(NULL), D1ij_(NULL)
{
    char vname[1024];

    strcpy(name_, "ft");

    nu_ = 1;

    int i = *carg;
    if (i >= narg)
        error->all(FLERR,"Missing crystal surface.");

    (*carg)++;
    crystal_surface_ = crystal_surface_factory(arg[i], narg, carg, arg, error);
    if (!crystal_surface_) {
        char errstr[1024];
        sprintf(errstr, "Could not find crystal surface '%s'.", arg[i]);
        error->all(FLERR,errstr);
    }

    if (*carg >= narg)
        error->all(FLERR,"Missing interaction model (force constants "
                   "definition).");

    i = *carg;
    (*carg)++;
    force_constants_ = force_constants_factory(arg[i], narg, carg, arg,
                                               crystal_surface_, force,
                                               memory, error);
    if (!force_constants_) {
        char errstr[1024];
        sprintf(errstr, "Could not find interaction model '%s'.", arg[i]);
        error->all(FLERR,errstr);
    }

    bool keyword_found = true;
    char *endptr;
    while (*carg < narg-1 && keyword_found) {
        keyword_found = false;
        if (!strcmp(arg[*carg], "height")) {
            (*carg)++;
            height_ = strtol(arg[*carg], &endptr, 10);
            if (endptr == arg[*carg]) {
                error_->all(FLERR,"Error converting height argument to "
                            "number.");
            }
            (*carg)++;
            keyword_found = true;
        }
    }

    ndof_ = 3*crystal_surface_->get_number_of_atoms();;
    dim_ = ndof_*nu_; // Tot dim=3*nc (in multiatom cell) or 3*num (1-atom cell)
}


FTStiffnessKernel::~FTStiffnessKernel()
{
    //delete [] linf_; // TAS3 does this prevent seg fault
    if (D0ii_)
        delete [] D0ii_;
    if (D1ii_)
        delete [] D1ii_;
    if (D0ij_)
        delete [] D0ij_;
    if (D1ij_)
        delete [] D1ij_;
}


/*!
 * Dump generic text info to file
 */
void FTStiffnessKernel::dump_info(FILE *f)
{
    StiffnessKernel::dump_info(f);

    crystal_surface_->dump_info(f);

    force_constants_->dump_info(f, crystal_surface_, error_);
}


/*!
 * Called before get_dynamical_matrices and used to tabulate Dij for each
 * pair ij.
 */
void FTStiffnessKernel::pre_compute()
{
    if (D0ii_)
        delete [] D0ii_;
    if (D1ii_)
        delete [] D1ii_;
    if (D0ij_)
        delete [] D0ij_;
    if (D1ij_)
        delete [] D1ij_;

    // Let the force constant module know about the crystal surface
    force_constants_->pre_compute(crystal_surface_, error_);

    // Total number of neighbors.
    int nn = crystal_surface_->get_number_of_neighbors();
    D0ij_ = new double[9*nn];
    D1ij_ = new double[9*nn];

    // Compute offsite elements.
    const CrystalSurface::lattice_neighbor_t *neigh = 
        crystal_surface_->get_neighbors();
    for (int i = 0; i < nn; i++, neigh++) {
        force_constants_->offsite(crystal_surface_, neigh, 0, &D0ij_[9*i],
                                  error_);
        force_constants_->offsite(crystal_surface_, neigh, 1, &D1ij_[9*i],
                                  error_);
    }

    // Number of atoms per surface unit cell.
    int nat = crystal_surface_->get_number_of_atoms();
    D0ii_ = new double[9*nat];
    D1ii_ = new double[9*nat];
    std::fill(D0ii_, D0ii_+9*nat, 0.0);
    std::fill(D1ii_, D1ii_+9*nat, 0.0);

    // Use sum rule to compute onsite elements.
    //   D_ii = \sum_j D_ij
    neigh = crystal_surface_->get_neighbors();
    double *D0ij = D0ij_;
    double *D1ij = D1ij_;
    for (int i = 0; i < nat; i++) {
#if 0
        mat<double> D0ii(3, &D0ii_[9*i]);
        mat<double> D1ii(3, &D1ii_[9*i]);

        int nn = crystal_surface_->get_number_of_neighbors(i);
        for (int j = 0; j < nn; j++, neigh++, D0ij+=9, D1ij+=9) {
            D0ii -= D0ij;
            D1ii -= D1ij;
        }

        mat<double,3> tmp0, tmp1;
        force_constants_->onsite(crystal_surface_, i, 0, tmp0.data(), error_);
        force_constants_->onsite(crystal_surface_, i, 1, tmp1.data(), error_);

        (D0ii - tmp0).print();
        (D1ii - tmp1).print();
#endif
        force_constants_->onsite(crystal_surface_, i, 0, &D0ii_[9*i], error_);
        force_constants_->onsite(crystal_surface_, i, 1, &D1ii_[9*i], error_);

#if GFMD_CHECK_POSITIVE_DEFINITE
        char name[1024];
        sprintf(name, "D0ii_[%i(nat=%i)]", i, nat);
        warn_positive_definite(name, 3, &D0ii_[9*i], error_);
        sprintf(name, "D1ii_[%i(nat=%i)]", i, nat);
        warn_positive_definite(name, 3, &D1ii_[9*i], error_);
#endif
    }
}


/*!
 * Function returns fac[] and D=[U0, ... V W] at this q vector  
 * Note: called each time the 'run' or 'minimize' command is called
 * Multi-atom unit cells written as single atom cells, per Appendix A.3
 */
void FTStiffnessKernel::get_dynamical_matrices(double qx, double qy,
                                               double_complex *xU0,
                                               double_complex *xU,
                                               double_complex *xV,
                                               double_complex dU)
{
    // Number of atoms per surface unit cell.
    int nat = crystal_surface_->get_number_of_atoms();

    assert(dim_ == 3*nat);

    mat<double_complex> U0(dim_, xU0);
    mat<double_complex> U(dim_, xU);
    mat<double_complex> V(dim_, xV);

    // Set all matrices to zero.
    U0.fill_with(0.0);
    U.fill_with(0.0);
    V.fill_with(0.0);

    mat<double, 3> cell(crystal_surface_->get_cell());

    // Add off-site terms to U matrices.
    int nn = crystal_surface_->get_number_of_neighbors();
    const CrystalSurface::lattice_neighbor_t *neigh =
        crystal_surface_->get_neighbors();
    for (int n = 0; n < nn; n++, neigh++) {
        // Get index of first atom.
        int i = neigh->indexi;

        // Get distance in number of cells. This determines the phase factor
        // below.
        double cx[3];
        crystal_surface_->cell_distance(i, *neigh, cx);

        // Compute phase.
        double_complex phase = exp(COMPLEX_NUMBER(0, cx[0]*qx+cx[1]*qy));

        // ab tells us whether it's the U(ab==0),V(ab==1),W,... matrix
        int ab = nearbyint(cx[2]);

        assert(std::abs(ab-cx[2]) < 1e-6);

        // Everything with ab < -1 and ab > 1 is ignored since those are
        // neighbors that are farther than one layer away. If the
        // interaction reaches this far the size of the unit cell should
        // be increased.
        if (ab >= -1 && ab <= 1) {
            // Get index of second atom.
            int j = neigh->indexj;

            // Get force constants.
            //force_constants_->offsite(crystal_surface_, neigh, Y.data());
            mat<double> Y(3, &D0ij_[9*n]);
            mat<double> Z(3, &D1ij_[9*n]);

            if (ab == 0) {
                // This needs to be added to the appropriate subblock of
                // U0 and U.

                for (int ii = 0; ii < 3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        U0[3*i+ii][3*j+jj] += phase*Y[ii][jj];
                        U[3*i+ii][3*j+jj] += phase*Z[ii][jj];
                    }
                }
            }
            else if (ab == -1) {
                // This needs to be added to the appropriate subblock of V.

                for (int ii = 0; ii < 3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        V[3*i+ii][3*j+jj] += phase*Z[ii][jj];
                    }
                }
            }
        }
    }

    // Add on-site terms to U0 and U.
    for (int i = 0; i < nat; i++) {
        //force_constants_->onsite(crystal_surface_, i, 0, Y.data());
        mat<double> Y(3, &D0ii_[9*i]);

        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                U0[3*i+ii][3*i+jj] += Y[ii][jj];
            }
        }

        //force_constants_->onsite(crystal_surface_, i, -1, Y.data());
        mat<double> Z(3, &D1ii_[9*i]);

        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                U[3*i+ii][3*i+jj] += Z[ii][jj];
            }
        }
    }

#ifdef GFMD_DEBUG
    if (std::abs(qx) < 1e-6 && std::abs(qy) < 1e-6) {
        double_complex trU0 = 0.0, trU = 0.0;
        for (int i = 0; i < dim_; i++) {
            trU0 += U0[i][i];
            trU += U[i][i];
        }
        std::cout << "USER-GFMD: DEBUG - trace[U0(q=0)] = " << creal(trU0)
                  << ", trace[U(q=0)] = " << creal(trU) << std::endl;
    }
#endif
}

void FTStiffnessKernel::get_force_at_gamma_point(double *f)
{
    force_constants_->linear(crystal_surface_, f, error_);

#if 0
    for (int i = 0; i < crystal_surface_->get_number_of_atoms(); i++) {
        const  double *r = crystal_surface_->get_position(i);
        std::cout << r[0] << ", " << r[1] << ", " << r[2] << ": " << f[i] << std::endl;
    }
#endif
}


} /* namespace */
