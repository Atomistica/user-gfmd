#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"
#include "math_const.h"

#include "fc_eam_fd.h"

FCEAMFD::FCEAMFD(int narg, int *carg, char **arg, CrystalSurface *surface, 
                 Force *force, Memory *memory, Error *error)
    : FCFiniteDifferences("eam/fd", surface, force, error, "eam/gf"),
      memory_(memory)
{
    manybody_ = true;

    pair_ = dynamic_cast<LAMMPS_NS::PairEAMGF*>(force->pair);
    pair_->init();

    surface->compute_neighbor_list(2*pair_->cutmax);
}


/*!
 * Dump generic text info to file
 */
void FCEAMFD::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    FCFiniteDifferences::dump_info(f, surf, error);
}


/*!
 * Compute forces on all atoms
 */
double FCEAMFD::energy_and_forces(const CrystalSurface *surf, int *mask, 
                                  double *u, double *f, Error *error)
{
  const CrystalSurface::lattice_neighbor_t *ii, *jj;
  int i,j,m,jnum,itype,jtype;
  double delx,dely,delz,esum,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;

  esum = 0.0;

  const double *x = surf->get_positions();
  const int *type = surf->get_types();
  int nlocal = surf->get_number_of_atoms();

  if (nlocal > pair_->nmax) {
    memory_->destroy(pair_->rho);
    memory_->destroy(pair_->fp);
    pair_->nmax = nlocal;
    memory_->create(pair_->rho,pair_->nmax,"pair:rho");
    memory_->create(pair_->fp,pair_->nmax,"pair:fp");
  }

  double *rho = pair_->rho;
  double *fp = pair_->fp;

  // zero out density and derivative of embedding energy

  for (i = 0; i < nlocal; i++) {
    rho[i] = 0.0;
    fp[i] = 0.0;
  }

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    jnum = surf->get_number_of_neighbors(i);

    jj = surf->get_neighbors(i);
    for (int njj = 0; njj < jnum; njj++, jj++) {
      j = jj->indexj;

      // pair eam uses half neighbor list
      if (i > j)  continue;

      delx = jj->r[0]+u[3*i+0]-u[3*j+0];
      dely = jj->r[1]+u[3*i+1]-u[3*j+1];
      delz = jj->r[2]+u[3*i+2]-u[3*j+2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < pair_->cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*pair_->rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,pair_->nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = pair_->rhor_spline[pair_->type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pair_->rhor_spline[pair_->type2rhor[itype][jtype]][m];
        rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      }
    }
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (i = 0; i < nlocal; i++) {
    // if mask is present, compute energy only for atoms with mask set
    if (mask && !mask[i])  continue;

    p = rho[i]*pair_->rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,pair_->nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = pair_->frho_spline[pair_->type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho[i] > pair_->rhomax) phi += fp[i] * (rho[i]-pair_->rhomax);
    esum += phi;
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    itype = type[i];

    delx = jj->r[0]+u[3*i+0]-u[3*j+0];
    dely = jj->r[1]+u[3*i+1]-u[3*j+1];
    delz = jj->r[2]+u[3*i+2]-u[3*j+2];
    jnum = surf->get_number_of_neighbors(i);

    jj = surf->get_neighbors(i);
    for (int njj = 0; njj < jnum; njj++, jj++) {
      j = jj->indexj;

      // pair eam uses half neighbor list
      if (i > j)  continue;

      delx = jj->r[0]+u[3*i+0]-u[3*j+0];
      dely = jj->r[1]+u[3*i+1]-u[3*j+1];
      delz = jj->r[2]+u[3*i+2]-u[3*j+2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < pair_->cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*pair_->rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,pair_->nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        coeff = pair_->rhor_spline[pair_->type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = pair_->rhor_spline[pair_->type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = pair_->z2r_spline[pair_->type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip; 
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip;

        double fac = 0;
        if (!mask)
          fac = 1.0;
        else if (mask[i] && mask[j])
          fac += 1.0;
        else if (mask[i] || mask[j])
          fac += 0.5;
        psip += fac*phip;
        esum += fac*phi;

#if 0
        if (!mask || (mask[i] && mask[j])) {
            psip += phip;
            esum += phi;
        }
#endif

        fpair = -psip*recip;

        f[3*i+0] += delx*fpair;
        f[3*i+1] += dely*fpair;
        f[3*i+2] += delz*fpair;
        f[3*j+0] -= delx*fpair;
        f[3*j+1] -= dely*fpair;
        f[3*j+2] -= delz*fpair;
      }
    }
  }

  return esum;
}
