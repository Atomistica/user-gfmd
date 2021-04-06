#ifdef FC_CLASS

FCStyle(eam,FCEAM)

#else

#ifndef __FC_EAM_H
#define __FC_EAM_H

#include "force_constants.h"

#include "pair_eam.h"

class FCEAM : public ForceConstants {
 public:
    FCEAM(int, int *, char **, CrystalSurface *, Force *, Error *);

    /*!
     * Off-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void offsite(const CrystalSurface *,
                         const CrystalSurface::lattice_neighbor_t *,
                         int, double *, Error *);

    /*!
     * On-site contribution to the force constant matrix, i.e. the
     * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
     */
    virtual void onsite(const CrystalSurface *, int, int, double *, Error *);

    /*!
     * More off-site contributions, from the Dij(A) [fourth term]
     */
    virtual void moreoffsite(const CrystalSurface *,
                             const CrystalSurface::lattice_neighbor_t *,
                             double *);

    /*!
     * First derivative of electron density function at
     * equilibrium lattice positions.
     */
    double dfdr(double, int, int);

    /*!
     * Linear force contribution for the free surface.
     */
    virtual void linear(const CrystalSurface *, double *, Error *error);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *, CrystalSurface *, Error *error);

    /*!
     * Set arrays of potential derivatives to be used in stiffness
     * calculations.  Arrays indexed by surface atom i and neighbor j.
     */
    void set_F_and_rho_from_paiream(const CrystalSurface *);

    /*!
     * First derivative of energy functional at equilibrium positions
     * Vector length is number of surface cell atoms.
     */
    double *equil_rho;
    double *dFdrho;
    double *d2Fdrho2;
    double **d2fdr2;
    double **equil_phip;
    double **equil_phipp;
    double **equil_kappa;

    /*!
     * Copies of the pair_eam variables
     */
     double cutmax;
     int nrho,nr;
     double **frho,**rhor,**z2r;
     int *type2frho,**type2rhor,**type2z2r;
     double dr,rdr,drho,rdrho,rhomax;
     double ***rhor_spline,***frho_spline,***z2r_spline;

 protected:
    /*!
     * The eam pair style class
     */
    LAMMPS_NS::PairEAM *pair_;
};

#endif

#endif
