#include <algorithm>
#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"
#include "math_const.h"

#include "fc_finite_differences.h"


FCFiniteDifferences::FCFiniteDifferences(const char *name,
                                         CrystalSurface *surface, 
                                         Force *force, Error *error,
                                         const char *pair_style)
    : ForceConstants(name, force, error, pair_style), delta_(1e-6),
      supercell_(NULL)
{
}


FCFiniteDifferences::~FCFiniteDifferences()
{
    if (supercell_)
        delete supercell_;
}


/*!
 * Dump generic text info to file
 */
void FCFiniteDifferences::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    ForceConstants::dump_info(f, surf, error);

    fprintf(f, "Using finite differences with a displacement of %f to compute "
            "the force-constant matrix.\n", delta_);
    fprintf(f, "Interaction cutoff = %f\n", get_cutoff());
}


/*!
 * Create a supercell for finite differences evaluation of the
 * force-constant matrix. Supercell needs to be at least 2*cutoff.
 */
CrystalSurface *FCFiniteDifferences::supercell(const CrystalSurface *surf,
                                               int &nx, int &ny, int &nz,
                                               Error *error)
{
    mat<double,3> cell(surf->get_cell());
    double lx, ly, lz;
    lx = cell[0][0];
    ly = cell[1][1];
    lz = sqrt(cell[0][2]*cell[0][2]+cell[1][2]*cell[1][2]+
              cell[2][2]*cell[2][2]);
    
    // create supercell of appropriate size (2*cutoff)
    CrystalSurface *supercell_ = new CrystalSurface(*surf);
    double cutoff = get_cutoff();
    nx = int(2*cutoff/lx)+1;
    ny = int(2*cutoff/ly)+1;
    nz = int(2*cutoff/lz)+1;
    if (nz <= 3) nz = 4;
    supercell_->reverse_indices();
    supercell_->supercell(nx, ny, nz, error);
    supercell_->reverse_indices();

#if 0
    for (int i = 0; i < surf->get_number_of_atoms(); i++) {
        std::cout << i << ": " << vec<double,3>(surf->get_position(i)) << std::endl;
    }
    for (int i = 0; i < supercell_->get_number_of_atoms(); i++) {
        std::cout << i << ": " << vec<double,3>(supercell_->get_position(i)) << std::endl;
    }
#endif

    return supercell_;
}


/*!
 * Determine which atoms will be handled by the GFMD part of the code.
 */
void FCFiniteDifferences::surface_mask(int nat, int nat0, int nz, int *mask)
{
    int nnz = nz*nat0;
    if (manybody_) {
        // for manybody potentials the energy for the upper half of the
        // cell is computed by the full potential
        for (int i = 0; i < nat; i++) {
            // this is the grid index, details depend on how supercell is
            // constructed in CrystalSurface::supercell
            int iu = i%nnz;
            mask[i] = iu < nat0/2 ? 0 : ( iu >= (nnz/2+nat0/2) ? 0 : 1 );
        }
    }
    else {
        for (int i = 0; i < nat; i++) {
            // this is the grid index, details depend on how supercell is
            // constructed in CrystalSurface::supercell
            int iu = i%nnz;
            mask[i] = iu >= nnz/2 ? 0 : 1;
        }
    }
}


/*!
 * Inform the force constants module about the lattice structure. Allow it
 * to construct supercell.
 */
void FCFiniteDifferences::pre_compute(const CrystalSurface *surf, Error *error)
{
    if (supercell_)
        delete supercell_;
    supercell_ = supercell(surf, nx_, ny_, nz_, error);   
}


/*!
 * Off-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
 */
void FCFiniteDifferences::offsite(const CrystalSurface *surf,
                                  const CrystalSurface::lattice_neighbor_t *pair,
                                  int ab, double *_D, Error *error)
{
    mat<double> D(3, _D);
    D = 0.0;

    // number of atoms in standard cell and supercell
    int nat0 = surf->get_number_of_atoms();
    int nat = supercell_->get_number_of_atoms();

    // mask for surface energies and forces
    int *mask = NULL;
#if 1
    if (!ab) {
        mask = new int[nat];
        surface_mask(nat, nat0, nz_, mask);
    }
#endif

    // displacements and forces
    double u[3*nat], f0[3*nat];
    double f1ik[3*nat], f2ik[3*nat], f1jk[3*nat], f2jk[3*nat];
    std::fill(u, u+3*nat, 0);

    // get indices in the supercell
    int i = pair->indexi;
    //std::cout << "nx = " << nx_ << ", ny = " << ny_ << ", nz = " << nz_ << ", modulo: " << modulo(pair->x, nx_) << ", " << modulo(pair->y, ny_) << ", " << modulo(pair->ab, nz_) << std::endl;
    int j = CrystalSurface::_idx(ny_, nz_, nat0, modulo(-pair->x, nx_),
                                 modulo(-pair->y, ny_), modulo(-pair->ab, nz_),
                                 pair->indexj);

    //std::cout << "offsite, ab = " << ab << " === " << i << "-" << j << " === " << *pair << std::endl;

#if 0
    vec<double,3> tmp = supercell_->get_position(i);
    tmp -= supercell_->get_position(j);
    std::cout << *pair << " - i = " << i << ", j = " << j << std::endl;
    std::cout << tmp << " " << vec<double,3>(pair->r) << " " << mat<double,3>(supercell_->get_cell()) << std::endl;
    std::cout << vec<double,3>(supercell_->get_position(i)) << " " << vec<double,3>(supercell_->get_position(j)) << " | " << vec<double,3>(surf->get_position(pair->indexi)) << " " << vec<double,3>(surf->get_position(pair->indexj)) << " | " << mat<double,3>(surf->get_cell()) << std::endl;
#endif

    std::fill(f0, f0+3*nat, 0);
    double e = energy_and_forces(supercell_, mask, u, f0, error);
    //printf("epot/at = %e\n", e);

    // wiggle in x, y and z direction
    for (int k = 0; k < 3; k++) {
        u[3*i+k] = -delta_;
        std::fill(f1ik, f1ik+3*nat, 0);
        double ei1 = energy_and_forces(supercell_, mask, u, f1ik, error);
        u[3*i+k] = delta_;
        std::fill(f2ik, f2ik+3*nat, 0);
        double ei2 = energy_and_forces(supercell_, mask, u, f2ik, error);
        u[3*i+k] = 0;

        u[3*j+k] = -delta_;
        std::fill(f1jk, f1jk+3*nat, 0);
        double ej1 = energy_and_forces(supercell_, mask, u, f1jk, error);
        u[3*j+k] = delta_;
        std::fill(f2jk, f2jk+3*nat, 0);
        double ej2 = energy_and_forces(supercell_, mask, u, f2jk, error);
        u[3*j+k] = 0;

#if 0
        if ((pair->indexi == 0 && pair->indexj == 2) ||
            (pair->indexi == 2 && pair->indexj == 0)) {
            std::cout << "xxx " << pair->indexi << "-" << pair->indexj << std::endl;
            for (int l = 0; l < 3*nat; l++) {
                vec<double,3> dr = supercell_->get_position(j);
                dr -= supercell_->get_position(l/3);
                std::cout << k << " " << (l%3) << " " << (l/3) << " " << dr.nrm2() << " " << dr << " " << fi1[l] << " " << fj1[l] << std::endl;
            }
        }
#endif

        for (int l = 0; l < 3; l++) {
            // finite differences approximation
            D[k][l] -= (f2ik[3*j+l]-f1ik[3*j+l])/(4*delta_);
            D[l][k] -= (f2jk[3*i+l]-f1jk[3*i+l])/(4*delta_);
#if 0
            if ((pair->indexi == 0 && pair->indexj == 2) ||
                (pair->indexi == 2 && pair->indexj == 0)) {
                printf("fj1: %e %e - %e %e\n", fi1[3*j+l], fi2[3*j+l], fj1[3*i+l], fj2[3*i+l]);
                printf("D[%i][%i] = %e\n", l, k, D[l][k]);
            }
#endif
        }

#ifdef GFMD_DEBUG
        double numfi = -(ei2-ei1)/(2*delta_);
        double numfj = -(ej2-ej1)/(2*delta_);
        //printf("%f %f %f %f %f\n", e, ei1, ei2, ej1, ej2);

        if (std::abs(numfi - f0[3*i+k]) > 1e-4) {
            printf("Warning: Analytical and numerical derivatives do not match "
                   "(atom %i, dir. %i) fi = %e, numfi = %e\n", i, k, f0[3*i+k],
                   numfi);        
        }
        if (std::abs(numfj - f0[3*j+k]) > 1e-4) {
            printf("Warning: Analytical and numerical derivatives do not match "
                   "(atom %i, dir. %i) fj = %e, numfj = %e\n", j, k, f0[3*j+k],
                   numfj);        
        }
#endif
    }

    if (mask)
        delete [] mask;

#if 0
    if (!ab) {
        std::cout << pair->indexi << " " << pair->indexj << std::endl;
        std::cout << D << std::endl;
    }
#endif
}


#if 0
/*!
 * On-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
 */
void FCFiniteDifferences::onsite(const CrystalSurface *surf, int i, int ab,
                                 double *_D, Error *error)
{
    mat<double> D(3, _D);
    D = 0.0;

    int n = surf->get_number_of_neighbors(i);
    const CrystalSurface::lattice_neighbor_t *neigh = surf->get_neighbors(i);

    for (int j = 0; j < n; j++, neigh++) {
        double tmp[9];
        offsite(surf, neigh, ab, tmp, error);
        D -= tmp;
    }
}
#endif


#if 1
/*!
 * On-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
 *
 * Note: ab==1 is bulk calculation, ab==0 surface
 */
void FCFiniteDifferences::onsite(const CrystalSurface *surf, int i, int ab,
                                 double *_D, Error *error)
{
    mat<double> D(3, _D);
    D = 0.0;

    // number of atoms in standard cell and supercell
    int nat0 = surf->get_number_of_atoms();
    int nat = supercell_->get_number_of_atoms();

    // mask for surface energies and forces
    int *mask = NULL;
    if (!ab) {
        mask = new int[nat];
        surface_mask(nat, nat0, nz_, mask);
    }

    // displacements and forces
    double u[3*nat], f0[3*nat];
    double f1[3*nat], f2[3*nat];
    std::fill(u, u+3*nat, 0);

    std::fill(f0, f0+3*nat, 0);
    double e = energy_and_forces(supercell_, mask, u, f0, error);

    // wiggle in x, y and z direction
    for (int k = 0; k < 3; k++) {
        u[3*i+k] = -delta_;
        std::fill(f1, f1+3*nat, 0);
        double ei1 = energy_and_forces(supercell_, mask, u, f1, error);
        u[3*i+k] = delta_;
        std::fill(f2, f2+3*nat, 0);
        double ei2 = energy_and_forces(supercell_, mask, u, f2, error);
        u[3*i+k] = 0;

        for (int l = 0; l < 3; l++) {
            // finite differences approximation
            D[k][l] -= (f2[3*i+l]-f1[3*i+l])/(4*delta_);
            D[l][k] -= (f2[3*i+l]-f1[3*i+l])/(4*delta_);
        }

#if 0
        if (i == 0) {
            printf("%i %i %i - %e %e - %e %e  %e %e  %e %e\n", i, k, ab, ei1-e, ei2-e,
                   f1[3*i+0], f2[3*i+0], f1[3*i+1], f2[3*i+1],
                   f1[3*i+2], f2[3*i+2]);
        }
#endif

#ifdef GFMD_DEBUG
        double numfi = -(ei2-ei1)/(2*delta_);
        //printf("%f %f %f %f %f\n", e, ei1, ei2, ej1, ej2);

        if (std::abs(numfi - f0[3*i+k]) > 1e-4) {
            printf("Warning: Analytical and numerical derivatives do not match "
                   "(atom %i, dir. %i) fi = %e, numfi = %e\n", i, k, f0[3*i+k],
                   numfi);        
        }
#endif
    }

    if (mask)
        delete [] mask;
}
#endif


/*!
 * Linear force contribution for the free surface.
 */
void FCFiniteDifferences::linear(const CrystalSurface *surf, double *linf,
                                 Error *error)
{
    int nx, ny, nz;
    CrystalSurface *supersurf = supercell(surf, nx, ny, nz, error);

    int nat0 = surf->get_number_of_atoms();
    int nat = supersurf->get_number_of_atoms();
    double u[3*nat], f[3*nat];
    int mask[nat];
    std::fill(u, u+3*nat, 0);
    std::fill(f, f+3*nat, 0);
    surface_mask(nat, nat0, nz, mask);

    double e = energy_and_forces(surf, NULL, u, f, error);

    std::cout << "USER-GFMD: Potential energy per atom in the ideal crystal = "
    	      << e/surf->get_number_of_atoms() << std::endl;
    double minval = *std::min_element(f, f+3*nat0);
    double maxval = *std::max_element(f, f+3*nat0);
    if (std::abs(minval) > 1e-9 || std::abs(maxval) > 1e-9) {
        std::cout << "USER-GFMD: Warning: Forces do not sum to zero in the "
                  << "bulk. (Min force = " << minval << ", max force = " 
                  << maxval << ".)" << std::endl;
    }

    std::fill(f, f+3*nat, 0);
    energy_and_forces(supersurf, mask, u, f, error);
#if 0
    for (int i = 0; i < nat; i++) {
        if (i < nat0) {
            std::cout << i << ": " << vec<double,3>(supersurf->get_position(i)) << " - " << mask[i] << " - " << vec<double,3>(&f[3*i]) << " | " << vec<double,3>(surf->get_position(i)) << std::endl;
        }
        else {
            std::cout << i << ": " << vec<double,3>(supersurf->get_position(i)) << " - " << mask[i] << " - " << vec<double,3>(&f[3*i]) << std::endl;
        }
    }
#endif

    for (int i = 0; i < nat0; i++) {
        linf[i] = f[3*i+2];
    }

#if 0
    double esum = 0.0;
    std::cout << "nat0 = " << nat0 << ", nz = " << nz << ", nat = " << nat << std::endl;
    for (int i = 0; i < nat; i++) {
        std::cout << i << " (" << mask[i] << "): " << f[3*i+2] << std::endl;
        esum += f[3*i+2];
    }
    std::cout << "esum = " << esum << std::endl;
#endif

    delete supersurf;
}
