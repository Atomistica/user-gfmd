#ifdef FC_CLASS

FCStyle(lj/cut,FCLJCut)

#else

#ifndef __FC_LJ_CUT_H
#define __FC_LJ_CUT_H

#include "force_constants.h"

#include "pair_lj_cut_gf.h"

class FCLJCut : public ForceConstants {
 public:
    FCLJCut(int, int *, char **, CrystalSurface *, Force *, Error *);

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
    virtual void linear(const CrystalSurface *, double *, Error *error);

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
     * The lj/cut pair style class
     */
    LAMMPS_NS::PairLJCutGF *pair_;
};

#endif

#endif
