#ifndef __FC_FINITE_DIFFERENCES_H
#define __FC_FINITE_DIFFERENCES_H

#include "force_constants.h"

class FCFiniteDifferences : public ForceConstants {
 public:
    FCFiniteDifferences(const char *, CrystalSurface *, Force *, Error *,
                        const char *pair_style=NULL);
    virtual ~FCFiniteDifferences();

    /*!
     * Inform the force constants module about the lattice structure. Allow it
     * to construct supercell.
     */
    virtual void pre_compute(const CrystalSurface *, Error *);

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *error);

    /*!
     * On-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
     */
    virtual void onsite(const CrystalSurface *, int, int, double *, Error *);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

    /*!
     * Create a supercell for finite differences evaluation of the
     * force-constant matrix. Supercell needs to be at least 2*cutoff.
     */
    virtual CrystalSurface *supercell(const CrystalSurface *, int &, int &,
                                      int &, Error *);

 protected:
    /*!
     * Return cutoff of potential
     */
    virtual double get_cutoff() = 0;

    /*!
     * Determine which atoms will be handled by the GFMD part of the code.
     */
    void surface_mask(int, int, int, int *);
    
    /*!
     * Compute energy and forces for a given atomic configuration
     */
    virtual double energy_and_forces(const CrystalSurface *, int *, double *,
                                     double *, Error *) = 0;

    /*!
     * Displacements for finite differences
     */
    double delta_;

    /*!
     * Supercell parameters
     */
    int nx_, ny_, nz_;

    /*!
     * The supercell
     */
    CrystalSurface *supercell_;
};

#endif
