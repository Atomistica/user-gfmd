#ifdef FC_CLASS

FCStyleWithMemory(eam/fd,FCEAMFD)

#else

#ifndef __FC_EAM_FD_H
#define __FC_EAM_FD_H

#include "fc_finite_differences.h"

#include "pair_eam_gf.h"

class FCEAMFD : public FCFiniteDifferences {
 public:
    FCEAMFD(int, int *, char **, CrystalSurface *, Force *, Memory *, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * Return cutoff of potential
     */
    virtual double get_cutoff() {
        return pair_->cutmax;
    }

    /*!
     * Compute energy and forces for a given atomic configuration
     */
    virtual double energy_and_forces(const CrystalSurface *, int *, double *,
                                     double *, Error *);

    /*!
     * The eam pair style class
     */
    LAMMPS_NS::PairEAMGF *pair_;

    /*!
     * Memory management
     */
    Memory *memory_;
};

#endif

#endif
