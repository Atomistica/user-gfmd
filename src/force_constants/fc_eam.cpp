#include <stdlib.h>
#include <string.h>
#include "linearalgebra.h"

#include "fc_eam.h"

FCEAM::FCEAM(int narg, int *carg, char **arg, CrystalSurface *surface,
             Force *force, Error *error)
    : ForceConstants("eam", force, error)
{
    pair_ = dynamic_cast<LAMMPS_NS::PairEAM*>(force->pair);
    pair_->init();

    // Copy pair EAM pointers to local copies
    // Allows one-to-one comparison of code blocks in fc_eam and pair_eam
    cutmax = pair_->cutmax;
    nrho = pair_->nrho;
    nr = pair_->nr;
    frho = pair_->frho;
    rhor = pair_->rhor;
    type2frho = pair_->type2frho;
    type2rhor = pair_->type2rhor;
    type2z2r = pair_->type2z2r;
    dr = pair_->dr;
    rdr = pair_->rdr;
    drho = pair_->drho;
    rdrho = pair_->rdrho;
    rhor_spline = pair_->rhor_spline;
    frho_spline = pair_->frho_spline;
    z2r_spline = pair_->z2r_spline;

    // Tristan, set cutoff here as xxx.
    surface->compute_neighbor_list(cutmax);
}


/*!
 * Dump generic text info to file
 */
void FCEAM::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    ForceConstants::dump_info(f, surf, error);
}


/*!
 * Off-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
 */
void FCEAM::offsite(const CrystalSurface *surf,
                    const CrystalSurface::lattice_neighbor_t *pair,
                    int ab, double *_D, Error *error)
{

    // No contribution from here if i and j are well-separated
    vec<double, 3> rji(pair->r);
    if (rji.nrm2() > cutmax) return;

    int i = pair->indexi; // surfcell_atom i
    int j = pair->indexj; // neighstruct_atom j
    int nat = surf->get_number_of_atoms();
    int nat_tlayer = nat/2;

    vec<double, 3> e(pair->r); // To atom j from i, eji
    e /= e.nrm2(); // Is this norm squared?

    mat<double> D(3, _D);
    outer(e, e, D);

    mat<double> perp = identity<double>(3) - D;

    // Embedded energy
    double k = dFdrho[i] * d2fdr2[i][j];
    double kappa = dFdrho[i] / (pair->rnorm) * dfdr(pair->rnorm,i,j);

    // k,kappa need twice these values if i and j are in the substrate 
    int ij_in_bs = 0;
    if (i > nat_tlayer-1) ij_in_bs++;
    if (i - pair->ab  > nat_tlayer-1) ij_in_bs++;

    k *= ij_in_bs;
    kappa *= ij_in_bs;

    // Pair potential EAM term.  carefully check (-) against
    // fcc100fteam_stiffness.cpp
    // Add only if i and j are in the b or s regions
    if (i > nat_tlayer-1) if (i - pair->ab > nat_tlayer-1) {
        k += equil_phipp[i][j];
        kappa += equil_kappa[i][j];
    }

    D *= -k;
    perp *= kappa;

    D -= perp;

    moreoffsite(surf, pair, _D);
}

/*!
 * Final term in the sum over j: atom j = atom i: the delta function term
 * "On-site" contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i != j.
 */
void FCEAM::onsite(const CrystalSurface *surf, int i, int ab, double *_D,
                   Error *error)
{
#if 0
    // No contribution from here if i and j are well-separated
    vec<double, 3> rji(pair->r);
    if (rji.nrm2() > cutmax) return;

    int i = pair->indexi; // surfcell_atom i
    int j = pair->indexj; // neighstruct_atom j

    vec<double, 3> e(pair->r); // To atom j from i
    e /= e.nrm2(); // Is this norm squared?

    mat<double> D(3, _D);
    outer(e, e, D);

    mat<double> perp = identity<double>(3) - D;

    // Embedded energy
    double k = dFdrho[i] * d2fdr2[i][j];
    double kappa = dFdrho[i] / pair->rnorm * dfdr(pair->rnorm,i,j);

    // D (always U in onsite) needs twice these values if j is in b or s
    int nat = surf->get_number_of_atoms();
    int nat_tlayer = nat/2;
    if (i - pair->ab  > nat_tlayer-1) {
      k *= 2.0;
      kappa *= 2.0;
    }

    // Pair potential.  Add if both atoms in substrate (b or s)
    if ((i > nat_tlayer-1) && (i - pair->ab > nat_tlayer-1)) {
        k += equil_phipp[i][j]; 
        kappa += equil_kappa[i][j];
    }

    D *= -k;
    perp *= kappa;

    D -= perp;
#endif
}

// First, Second, Third terms cancel by symmetry of the FCC100 lattice.
// Copy those over once FCC100 checks out.

// Fourth Dij_A term
void FCEAM::moreoffsite(const CrystalSurface *surf,
                    const CrystalSurface::lattice_neighbor_t *pair,
                    double *_D)
{
    double dfdrki, dfdrkj;
    int ktype;

    int nat = surf->get_number_of_atoms();
    int nat_tlayer = nat/2;

    int i = pair->indexi; // a surface cell atom
    int j = pair->indexj; // a neighbor struct atom

    vec<double, 3> eki;
    vec<double, 3> ekj;
    mat<double> eki_ekj  = identity<double>(3);
    mat<double> D(3, _D);

    int itype = surf->get_type(i);
    int jtype = pair->jtype;

    // To atom j from i, rji
    vec<double, 3> rji(pair->r);

    const CrystalSurface::lattice_neighbor_t *neighk = 
        surf->get_neighbors(i);

    // Loop over neighbors of i and check if within rcut of i,j
    for (int k = 0; k < surf->get_number_of_neighbors(i); k++) {

        vec<double, 3> rki(neighk->r);
        vec<double, 3> rkj = rki - rji;

        // Only include if k is in the substrate.
        // Derivatives are zero if rki or rkj > cutmax.
        if ((rki.nrm2() <= cutmax) && (rkj.nrm2() <= cutmax)) {

            eki = rki;
            eki /= eki.nrm2();
            ekj = rkj;
            ekj /= ekj.nrm2();
            outer(eki, ekj, eki_ekj);

            if (i - neighk->ab  > nat_tlayer-1) {

                ktype = neighk->jtype;

                dfdrki = dfdr(rki.nrm2(),itype,ktype); // [i,k];
                dfdrkj = dfdr(rkj.nrm2(),jtype,ktype); // Check signs

                D += eki_ekj *dfdrki * dfdrkj * d2Fdrho2[i]; // watch casting
            }
        }
        neighk++;
    }
}

/*!
 * For each atom i in the surface cell,
 * use pair_eam.cpp splines to reference the F(rho) function.
 * Positions are given for the ideal lattice from CrystalSurface.
 */
void FCEAM::set_F_and_rho_from_paiream(const CrystalSurface *surf)
{
    double p,phi,z2,z2p,z2pp,equil_r,recip;
    double *coeff;
    int m, maxneighs = 0;
    int nat = surf->get_number_of_atoms();

    for (int i=0; i<nat; i++) 
        maxneighs = MAX(maxneighs,surf->get_number_of_neighbors(i));

    // Allocate memory for the 1D and 2D arrays of doubles
    equil_rho  = new double[nat];
    dFdrho     = new double[nat];
    d2Fdrho2   = new double[nat];
    d2fdr2     = new double*[nat];
    equil_phip = new double*[nat];
    equil_phipp= new double*[nat];
    equil_kappa= new double*[nat];
    for (int i = 0; i < nat; i++) {
        d2fdr2[i]     = new double[maxneighs];
        equil_phip[i] = new double[maxneighs];
        equil_phipp[i]= new double[maxneighs];
        equil_kappa[i]= new double[maxneighs];
    }
    // Should check if any are NULL and delete again
    //mm->cr(equil_rho,nat,"FCEAM:equil_rho");
    //    equil_rho = (double *) malloc(nat*sizeof(double));
    //mm->cr(dFdrho,nat,"FCEAM:dFdrho");
    //    dFdrho = (double *) malloc(nat*sizeof(double));
    //mm->cr(d2Fdrho2,nat,"FCEAM:d2Fdrho2");
    //    d2Fdrho2 = (double *) malloc(nat*sizeof(double));
    //mm->cr(d2fdr2,nat,maxneighs,"FCEAM:d2fdr2");
    //    d2fdr2 = (double *) malloc(nat*sizeof(double));
    //mm->cr(equil_phip,nat,maxneighs,"FCEAM:equil_phip");
    //    equil_phip = (double *) malloc(nat*sizeof(double));
    //mm->cr(equil_phipp,nat,maxneighs,"FCEAM:equil_phipp");
    //    equil_phipp = (double *) malloc(nat*sizeof(double));
    //mm->cr(equil_kappa,nat,maxneighs,"FCEAM:equil_kappa");
    //    equil_kappa = (double *) malloc(nat*sizeof(double));

    std::fill(equil_rho, equil_rho+nat, 0.0);
    std::fill(dFdrho, dFdrho+nat, 0.0);
    std::fill(d2Fdrho2, d2Fdrho2+nat, 0.0);
    for (int i = 0; i < nat; i++) {
        for (int j = 0; j < surf->get_number_of_neighbors(i); j++) {
            equil_phip[i][j]  = 0.0;
            equil_phipp[i][j] = 0.0;
            equil_kappa[i][j] = 0.0;
            d2fdr2[i][j] = 0.0; // equil_rhojpp
        }
    }

    // Similar to pair_eam compute() code.
    // But instead of looping over atoms, loop over latticeneighbor structs.

    // Fill in equil_rho[i] and F(rho) derivatives at neighbor distances.  
    // Then fill in equil_fp[i].

    // Adapted from the loop at line 169 in pair_eam.cpp
    // Loop each surface cell atom i
    for (int i = 0; i < nat; i++) {
        int itype = surf->get_type(i);

	// Loop through neighbors j of surface cell atom i
        const CrystalSurface::lattice_neighbor_t *neighj = 
            surf->get_neighbors(i);
        for (int j = 0; j < surf->get_number_of_neighbors(i); j++) {

            // Set values from all array that are functions of distance
            equil_r = neighj->rnorm; 
            if (equil_r < cutmax) {
                int jtype = neighj->jtype;

                // p - Interpolation value for distance between i and j
                p = (equil_r)*rdr + 1.0;
                m = static_cast<int> (p);
                m = MIN(m,nr-1);
                p -= m;
                p = MIN(p,1.0);

                // Get pair potential derivative
                coeff = z2r_spline[type2z2r[itype][jtype]][m];
                z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
                z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
                recip = 1.0/equil_r;
                phi = z2*recip;
                // Finally, the derivative
                equil_phip[i][j] = z2p*recip - phi*recip;

                // And second deriv is then: z2pp/r - z2p/r/r - phip/r + phi/r/r
                z2pp = 2.0*coeff[0]*p + coeff[1];
                equil_phipp[i][j] = recip*(z2pp - z2p*recip - equil_phip[i][j] +
                                           phi*recip);

                // Force constant for bond rotation at equilibrium
                // kappa = kperp = -f/r0 = dV(r0)/dr / r0
                equil_kappa[i][j] = equil_phip[i][j] * recip;

                // contribution to density at atom i from atom j
                coeff = rhor_spline[type2rhor[jtype][itype]][m];
                equil_rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + 
                    coeff[6];

                // // density derivative at atom i from atom j
                // // (If index by jshell instead of j, just overwrite)
                // coeff = rhor_spline[type2rhor[jtype][itype]][m];
                // equil_rhojp[i][j] = (coeff[0]*p + coeff[1])*p + coeff[2];

                // density 2nd derivative at atom i from atom j
                d2fdr2[i][j] = 2.0*coeff[0]*p + coeff[1]; // equil_rhojpp

            }
            neighj++;
        }
    }
    
    // This mimics the loop at line 213 in pair_eam.cpp
    // Interpolate for derivative of embedding energy(rho) at each surface
    // atom i: fp[i]
    for (int i = 0; i < nat; i++) {
        int itype = surf->get_type(i);
        p = equil_rho[i]*rdrho + 1.0;
        m = static_cast<int> (p);
        m = MAX(1,MIN(m,nrho-1));  
        p -= m;
        p = MIN(p,1.0);
        coeff = frho_spline[type2frho[itype]][m];
        dFdrho[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
        d2Fdrho2[i] = 2.0*coeff[0]*p + coeff[1];
    }
}

/*!
 * Derivative of the electron density function at given distance
 */
double FCEAM::dfdr(double distr, int itype, int jtype)
{
    if (distr > cutmax) return 0;

    int m;
    double p;
    double *coeff;

    p = sqrt(distr)*rdr + 1.0;
    m = static_cast<int> (p);
    m = MIN(m,nr-1);
    p -= m;
    p = MIN(p,1.0);

    coeff = rhor_spline[type2rhor[jtype][itype]][m];
    return (coeff[0]*p + coeff[1])*p + coeff[2];
}

/*!
 * Linear force contribution for the free surface.
 */
void FCEAM::linear(const CrystalSurface *surf, double *linf, Error *error)
{
    double zfraction, embedforce, pairforce;

    // number of atoms in the surface cell (t and b layers)
    int nat = surf->get_number_of_atoms();

    // number of atoms in the surface cell t layer
    int nat_tlayer = nat/2;

    std::fill(linf, linf+nat, 0.0);

    // Loop through atoms of surface cell
    // Assume i < nat_tlayer are within the t layer
    // (Cannot arbitrarily assign order of atoms within the surface cell: 
    //   ab and i must refer to the same structure)
    for (int i = 0; i < nat; i++) {
        int itype = surf->get_type(i);

        const CrystalSurface::lattice_neighbor_t *neighj = 
            surf->get_neighbors(i);

	// Loop through neighbors j of surface cell atom i
        for (int j = 0; j < surf->get_number_of_neighbors(i); j++) {

	    // fraction of vector in z-direction: -eij * ez
	    zfraction = (neighj->r)[2] / neighj->rnorm;

	    // Embedded contribution to linf[i] is d/dui of Eb+Es at equilibrium.
	    embedforce = dFdrho[i] * dfdr(neighj->rnorm,i,j) * zfraction; // dfdr=equil_rhojp
	    if (i > nat_tlayer-1) linf[i] += embedforce;
	    if (i-neighj->ab > nat_tlayer-1) linf[i] += embedforce;

	    // Pair potential contribution to linf[i] is d/dui of Ebs at equilibrium.
	    pairforce = equil_phip[i][j] * zfraction; // dphidr=equil_phip
	    if ((i > nat_tlayer-1) && (i-neighj->ab > nat_tlayer-1)) 
	        linf[i] += pairforce; // Not the negative TAS?

            neighj++;
        }
    }

    // Subtract net z-force
    double linfsum = 0;
    for (int i = 0; i < nat; i++) linfsum  += linf[i];
    for (int i = 0; i < nat; i++) linf[i] -= (linfsum / nat);
}

