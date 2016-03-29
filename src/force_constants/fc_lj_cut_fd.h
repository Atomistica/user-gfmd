#ifdef FC_CLASS

FCStyle(lj/cut/fd,FCLJCutFD)

#else

#ifndef __FC_LJ_CUT_FD_H
#define __FC_LJ_CUT_FD_H

#include "fc_finite_differences.h"

#ifdef GFMD_DEBUG
#include "fc_lj_cut.h"
#endif

#include "pair_lj_cut_gf.h"

class FCLJCutFD : public FCFiniteDifferences {
 public:
    FCLJCutFD(int, int *, char **, CrystalSurface *, Force *, Error *);
#ifdef GFMD_DEBUG
    virtual ~FCLJCutFD();

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *);

    /*!
     * On-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
     */
    virtual void onsite(const CrystalSurface *, int, int, double *, Error *);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *, Error *);
#endif

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * Return cutoff of potential
     */
    virtual double get_cutoff() {
        return pair_->cut_global;
    }

    /*!
     * Compute energy and forces for a given atomic configuration
     */
    virtual double energy_and_forces(const CrystalSurface *, int *, double *,
                                     double *, Error *);

    /*!
     * The lj/cut pair style class
     */
    LAMMPS_NS::PairLJCutGF *pair_;

#ifdef GFMD_DEBUG
    FCLJCut *debug_fc_lj_cut_;
#endif
};

#endif

#endif
