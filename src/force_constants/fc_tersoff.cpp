#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"
#include "math_const.h"

#include "fc_tersoff.h"

FCTersoff::FCTersoff(int narg, int *carg, char **arg, CrystalSurface *surface, 
                     Force *force, Error *error)
    : FCFiniteDifferences("tersoff", surface, force, error, "tersoff/gf")
{
    manybody_ = true;

    pair_ = dynamic_cast<LAMMPS_NS::PairTersoffGF*>(force->pair);
    pair_->init();

    surface->compute_neighbor_list(get_cutoff());
}


/*!
 * Dump generic text info to file
 */
void FCTersoff::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    FCFiniteDifferences::dump_info(f, surf, error);
}


/*!
 * Compute forces on all atoms
 */
double FCTersoff::energy_and_forces(const CrystalSurface *surf, int *mask, 
                                    double *u, double *f, Error *error)
{
    const CrystalSurface::lattice_neighbor_t *ii, *jj, *kk;
    int i,j,k,jnum;
    int itype,jtype,ktype,iparam_ij,iparam_ijk;
    double delx,dely,delz,esum,evdwl,fpair;
    double rsq,rsq1,rsq2;
    double delr1[3],delr2[3],fi[3],fj[3],fk[3];
    double zeta_ij,prefactor;
    
    esum = 0.0;

    const double *x = surf->get_positions();
    const int *type = surf->get_types();
    int nlocal = surf->get_number_of_atoms();

    // loop over full neighbor list of my atoms

    for (i = 0; i < nlocal; i++) {
        // if mask is present, compute only for atoms with mask set
        if (mask && !mask[i])  continue;

        itype = pair_->map[type[i]];

        // two-body interactions

        jnum = surf->get_number_of_neighbors(i);
        jj = surf->get_neighbors(i);
        for (int njj = 0; njj < jnum; njj++, jj++) {
            j = jj->indexj;

            // if mask is present, compute only for atoms with mask set for
            // both atoms
            if (mask && !mask[j])  continue;

            jtype = pair_->map[type[j]];

            delx = jj->r[0]+u[3*i+0]-u[3*j+0];
            dely = jj->r[1]+u[3*i+1]-u[3*j+1];
            delz = jj->r[2]+u[3*i+2]-u[3*j+2];
            rsq = delx*delx + dely*dely + delz*delz;

            iparam_ij = pair_->elem3param[itype][jtype][jtype];
            if (rsq > pair_->params[iparam_ij].cutsq) continue;

            pair_->repulsive(&pair_->params[iparam_ij],rsq,fpair,1,evdwl);

            // double counting term
            evdwl *= 0.5;
            fpair *= 0.5;

            f[3*i+0] += delx*fpair;
            f[3*i+1] += dely*fpair;
            f[3*i+2] += delz*fpair;
            f[3*j+0] -= delx*fpair;
            f[3*j+1] -= dely*fpair;
            f[3*j+2] -= delz*fpair;

            esum += evdwl;
        }

        // three-body interactions
        // skip immediately if I-J is not within cutoff

        jj = surf->get_neighbors(i);
        for (int njj = 0; njj < jnum; njj++, jj++) {
            j = jj->indexj;
            jtype = pair_->map[type[j]];
            iparam_ij = pair_->elem3param[itype][jtype][jtype];

            delr1[0] = -jj->r[0]-u[3*i+0]+u[3*j+0];
            delr1[1] = -jj->r[1]-u[3*i+1]+u[3*j+1];
            delr1[2] = -jj->r[2]-u[3*i+2]+u[3*j+2];
            rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
            if (rsq1 > pair_->params[iparam_ij].cutsq) continue;

            // accumulate bondorder zeta for each i-j interaction via loop
            // over k

            zeta_ij = 0.0;

            kk = surf->get_neighbors(i);
            for (int nkk = 0; nkk < jnum; nkk++, kk++) {
                if (jj == kk) continue;
                k = kk->indexj;
                ktype = pair_->map[type[k]];
                iparam_ijk = pair_->elem3param[itype][jtype][ktype];
                
                delr2[0] = -kk->r[0]-u[3*i+0]+u[3*k+0];
                delr2[1] = -kk->r[1]-u[3*i+1]+u[3*k+1];
                delr2[2] = -kk->r[2]-u[3*i+2]+u[3*k+2];
                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + 
                    delr2[2]*delr2[2];
                if (rsq2 > pair_->params[iparam_ijk].cutsq) continue;

                zeta_ij += pair_->zeta(&pair_->params[iparam_ijk],rsq1,rsq2,
                                       delr1,delr2);
            }

            // pairwise force due to zeta

            pair_->force_zeta(&pair_->params[iparam_ij],rsq1,zeta_ij,fpair,
                              prefactor,1,evdwl);

            f[3*i+0] += delr1[0]*fpair;
            f[3*i+1] += delr1[1]*fpair;
            f[3*i+2] += delr1[2]*fpair;
            f[3*j+0] -= delr1[0]*fpair;
            f[3*j+1] -= delr1[1]*fpair;
            f[3*j+2] -= delr1[2]*fpair;

            esum += evdwl;

            // attractive term via loop over k

            kk = surf->get_neighbors(i);
            for (int nkk = 0; nkk < jnum; nkk++, kk++) {
                if (jj == kk) continue;
                k = kk->indexj;
                ktype = pair_->map[type[k]];
                iparam_ijk = pair_->elem3param[itype][jtype][ktype];

                delr2[0] = -kk->r[0]-u[3*i+0]+u[3*k+0];
                delr2[1] = -kk->r[1]-u[3*i+1]+u[3*k+1];
                delr2[2] = -kk->r[2]-u[3*i+2]+u[3*k+2];
                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + 
                    delr2[2]*delr2[2];
                if (rsq2 > pair_->params[iparam_ijk].cutsq) continue;

                pair_->attractive(&pair_->params[iparam_ijk],prefactor,
                                  rsq1,rsq2,delr1,delr2,fi,fj,fk);

                f[3*i+0] += fi[0];
                f[3*i+1] += fi[1];
                f[3*i+2] += fi[2];
                f[3*j+0] += fj[0];
                f[3*j+1] += fj[1];
                f[3*j+2] += fj[2];
                f[3*k+0] += fk[0];
                f[3*k+1] += fk[1];
                f[3*k+2] += fk[2];
            }
        }
    }

    return esum; 
}
