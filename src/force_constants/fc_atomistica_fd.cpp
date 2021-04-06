#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"
#include "math_const.h"

#include "fc_atomistica_fd.h"

FCAtomisticaFD::FCAtomisticaFD(int narg, int *carg, char **arg,
                               CrystalSurface *surface, Force *force,
                               Error *error)
    : FCFiniteDifferences("atomistica/fd", surface, force, error, "atomistica"),
      nx_(-1), ny_(-1), nz_(-1), supercell_(NULL), nmax_(-1), nneighbmax_(-1),
      r_(NULL), f_(NULL), m_(NULL), tag_(NULL), neighb_(NULL), seed_(NULL),
      last_(NULL)
{
    manybody_ = true;

    pair_ = dynamic_cast<LAMMPS_NS::PairAtomistica*>(force->pair);
    pair_->init();

    surface->compute_neighbor_list(get_cutoff());
}


FCAtomisticaFD::~FCAtomisticaFD()
{
    if (supercell_) delete supercell_;

    if (r_) delete [] r_;
    if (f_) delete [] r_;
    if (m_) delete [] r_;
    if (tag_) delete [] r_;
    if (seed_) delete [] r_;
    if (last_) delete [] r_;
    if (neighb_) delete [] neighb_;
}


/*!
 * Dump generic text info to file
 */
void FCAtomisticaFD::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    FCFiniteDifferences::dump_info(f, surf, error);
}


/*!
 * Compute forces on all atoms
 */
double FCAtomisticaFD::energy_and_forces(const CrystalSurface *surf, int *mask, 
                                         double *u, double *f, Error *error)
{
    if (supercell_) {
        if (supercell_->get_number_of_atoms()/(nx_*ny_*nz_) != 
                surf->get_number_of_atoms() ||
            supercell_->get_number_of_neighbors()/(nx_*ny_*nz_) != 
                surf->get_number_of_neighbors()) {
            delete supercell_;
            supercell_ = NULL;
        }
    }

    if (!supercell_) {
        const double *cell = surf->get_cell();
        double cut = get_cutoff();
        int cx, cy, cz;
        cx = int(cut/cell[0])+1;
        cy = int(cut/cell[4])+1;
        cz = int(cut/cell[8])+1;
        nx_ = 2*cx+1;
        ny_ = 2*cy+1;
        nz_ = 2*cz+1;

        supercell_ = new CrystalSurface(*surf);
        supercell_->supercell(nx_, ny_, nz_, error, cx, cy, cz);
    }

    const double *x = supercell_->get_positions();
    int *type = supercell_->get_types();
    int nlocal = surf->get_number_of_atoms();
    int nall = supercell_->get_number_of_atoms();
    int nneighb = supercell_->get_number_of_neighbors();

    if (nmax_ < nall) {
        if (r_) delete [] r_;
        if (f_) delete [] r_;
        if (m_) delete [] r_;
        if (tag_) delete [] r_;
        if (seed_) delete [] r_;
        if (last_) delete [] r_;

        r_ = new double[3*nall];
        f_ = new double[3*nall];
        m_ = new int[nall];
        tag_ = new int[nall];
        seed_ = new intptr_t[nall];
        last_ = new intptr_t[nall];

        nmax_ = nall;
    }

    if (nneighbmax_ < nneighb) {
        if (neighb_) delete [] neighb_;

        neighb_ = new int[nneighb];

        nneighbmax_ = nneighb;
    }

    // ccopy particle positions and convert neighbor list
    nneighb = 0;
    for (int i = 0; i < nall; i++) {
        int tagi = i%nlocal;
        if (mask) {
            m_[i] = mask[tagi];
        }
        else {
            m_[i] = 1;
        }
        tag_[i] = tagi;
        r_[3*i+0] = x[3*i+0] + u[3*tagi+0];
        r_[3*i+1] = x[3*i+1] + u[3*tagi+1];
        r_[3*i+2] = x[3*i+2] + u[3*tagi+2];

        seed_[i] = nneighb+1;

        int jnum = supercell_->get_number_of_neighbors(i);
        const CrystalSurface::lattice_neighbor_t *jj = 
            supercell_->get_neighbors(i);
        for (int njj = 0; njj < jnum; njj++, jj++, nneighb++) {
            neighb_[nneighb] = jj->indexj;
        }

        last_[i] = nneighb;
    }

    assert(nneighb == supercell_->get_number_of_neighbors());

    // set pointers in particles object
    particles_set_pointers(pair_->particles_,nall,nlocal,nall,tag_,type,r_);

    // set pointers in neighbor list object
    neighbors_set_pointers(pair_->neighbors_,nall,seed_,last_,nneighb,neighb_);

    // compute energy and forces
    double epot = 0.0, wpot[9];
    std::fill(wpot, wpot+9, 0.0);
    std::fill(f_, f_+3*nall, 0.0);

    int ierror;
    pair_->class_->energy_and_forces(pair_->potential_,pair_->particles_,
                                     pair_->neighbors_,&epot,f_,wpot,m_,
                                     NULL,NULL,&ierror);
    error2lmp(error,FLERR,ierror);

    // sum forces to output buffer
    for (int i = 0; i < nall; i++) {
        int tagi = i%nlocal;
        f[3*tagi+0] += f_[3*i+0];
        f[3*tagi+1] += f_[3*i+1];
        f[3*tagi+2] += f_[3*i+2];
    }

    return epot;
}
