#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"

#include "fc_lj_smooth.h"

FCLJSmooth::FCLJSmooth(int narg, int *carg, char **arg, CrystalSurface *surface,
                       Force *force, Error *error)
    : ForceConstants("lj/smooth", force, error)
{
    pair_ = dynamic_cast<LAMMPS_NS::PairLJSmooth*>(force->pair);
    pair_->init();

    surface->compute_neighbor_list(pair_->cut_global);
}


/*!
 * Dump generic text info to file
 */
void FCLJSmooth::dump_info(FILE *f, CrystalSurface *surf, Error *error)
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
void FCLJSmooth::offsite(const CrystalSurface *surf,
                         const CrystalSurface::lattice_neighbor_t *pair,
                         int ab, double *_D, Error *)
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
double FCLJSmooth::first_derivative(int itype, int jtype, double r)
{
    double rsq = r*r;
    double r2inv = 1.0/rsq;
    if (rsq < pair_->cut_inner_sq[itype][jtype]) {
        double r6inv = r2inv*r2inv*r2inv;
        return r6inv * (pair_->lj2[itype][jtype] -
                        pair_->lj1[itype][jtype]*r6inv) / r;
    } else if (rsq < pair_->cut[itype][jtype]*pair_->cut[itype][jtype]) {
        double t = r - pair_->cut_inner[itype][jtype];
        double tsq = t*t;
        return -(pair_->ljsw1[itype][jtype] + 
                 pair_->ljsw2[itype][jtype]*t + 
                 pair_->ljsw3[itype][jtype]*tsq + 
                 pair_->ljsw4[itype][jtype]*tsq*t);
    }
    return 0.0;
}


/*!
 * Second derivative of pair expression
 */
double FCLJSmooth::second_derivative(int itype, int jtype, double r)
{
    double rsq = r*r;
    double r2inv = 1.0/rsq;
    if (rsq < pair_->cut_inner_sq[itype][jtype]) {
        double r6inv = r2inv*r2inv*r2inv;
        return r6inv * r2inv * (13*pair_->lj1[itype][jtype]*r6inv-
                                7*pair_->lj2[itype][jtype]);
    } else if (rsq < pair_->cut[itype][jtype]*pair_->cut[itype][jtype]) {
        double t = r - pair_->cut_inner[itype][jtype];
        double tsq = t*t;
        return -(pair_->ljsw2[itype][jtype] + 
                 2*pair_->ljsw3[itype][jtype]*t +
                 3*pair_->ljsw4[itype][jtype]*tsq);
    }
    return 0.0;
}


/*!
 * Linear force contribution for the free surface.
 */
void FCLJSmooth::linear(const CrystalSurface *surf, double *linf, Error *error)
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
            if (neigh->ab - i <= 0) { // TAS added - i
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
