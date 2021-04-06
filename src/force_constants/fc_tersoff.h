#ifdef FC_CLASS

FCStyle(tersoff,FCTersoff)

#else

#ifndef __FC_TERSOFF_H
#define __FC_TERSOFF_H

#include "fc_finite_differences.h"

#include "pair_tersoff_gf.h"

class FCTersoff : public FCFiniteDifferences {
 public:
    FCTersoff(int, int *, char **, CrystalSurface *, Force *, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * Return cutoff of potential
     */
    virtual double get_cutoff() {
        return 2*pair_->cutmax;
    }

    /*!
     * Compute energy and forces for a given atomic configuration
     */
    virtual double energy_and_forces(const CrystalSurface *, int *, double *,
                                     double *, Error *);

    /*!
     * The eam pair style class
     */
    LAMMPS_NS::PairTersoffGF *pair_;
};

#endif

#endif
