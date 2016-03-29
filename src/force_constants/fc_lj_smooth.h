#ifdef FC_CLASS

FCStyle(lj/smooth,FCLJSmooth)

#else

#ifndef __FC_LJ_SMOOTH_H
#define __FC_LJ_SMOOTH_H

#include "force_constants.h"

#include "pair_lj_smooth.h"

class FCLJSmooth : public ForceConstants {
 public:
    FCLJSmooth(int, int *, char **, CrystalSurface *, Force *, Error *);

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * First derivative of pair expression
     */
    virtual double first_derivative(int, int, double);

    /*!
     * Second derivative of pair expression
     */
    virtual double second_derivative(int, int, double);

    /*!
     * The lj/smooth pair style class
     */
    LAMMPS_NS::PairLJSmooth *pair_;
};

#endif

#endif
