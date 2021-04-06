#ifdef FC_CLASS

FCStyle(atomistica/fd,FCAtomisticaFD)

#else

#ifndef __FC_ATOMISTICA_FD_H
#define __FC_ATOMISTICA_FD_H

#include "fc_finite_differences.h"

#include "pair_atomistica.h"

class FCAtomisticaFD : public FCFiniteDifferences {
 public:
    FCAtomisticaFD(int, int *, char **, CrystalSurface *, Force *, Error *);
    virtual ~FCAtomisticaFD();

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *);

 protected:
    /*!
     * Return cutoff of potential
     */
    virtual double get_cutoff() {
        return pair_->rangemax_;
    }

    /*!
     * Compute energy and forces for a given atomic configuration
     */
    virtual double energy_and_forces(const CrystalSurface *, int *, double *,
                                     double *, Error *);

    /*!
     * The eam pair style class
     */
    LAMMPS_NS::PairAtomistica *pair_;

    /*!
     * Number of supercells
     */
    int nx_, ny_, nz_;

    /*!
     * Supercell including ghost particle region
     */
    CrystalSurface *supercell_;

    /*!
     * Length of buffers
     */
    int nmax_, nneighbmax_;

    /*!
     * Temporary buffers
     */
    double *r_, *f_;
    int *m_, *tag_, *neighb_;
    intptr_t *seed_, *last_;

    /*!
     * Memory management
     */
    Memory *memory_;
};

#endif

#endif
