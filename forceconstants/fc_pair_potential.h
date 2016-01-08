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
