/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <iostream>

#include "pair_tersoff_gf.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "suffix.h"
#include "tokenizer.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

PairTersoffGF::PairTersoffGF(LAMMPS *lmp) : PairTersoff(lmp)
{
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoffGF::~PairTersoffGF()
{
}

/* ---------------------------------------------------------------------- */

void PairTersoffGF::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  if (shift_flag) {
    if (evflag) {
      if (eflag) {
        if (vflag_atom) eval<1,1,1,1>();
        else eval<1,1,1,0>();
      } else {
        if (vflag_atom) eval<1,1,0,1>();
        else eval<1,1,0,0>();
      }
    } else eval<1,0,0,0>();

  } else {

    if (evflag) {
      if (eflag) {
        if (vflag_atom) eval<0,1,1,1>();
        else eval<0,1,1,0>();
      } else {
        if (vflag_atom) eval<0,1,0,1>();
        else eval<0,1,0,0>();
      }
    } else eval<0,0,0,0>();
  }
}

template <int SHIFT_FLAG, int EVFLAG, int EFLAG, int VFLAG_ATOM>
void PairTersoffGF::eval()
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double fforce;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double r1_hat[3],r2_hat[3];
  double zeta_ij,prefactor;
  double forceshiftfac;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *gflag = atom->gflag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

  if (!atom->gfmd_flag) gflag = nullptr;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      // if mask is present in both, do not compute those interactions
      if (gflag && gflag[i] && gflag[j]) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // shift rsq and store correction for force

      if (SHIFT_FLAG) {
        double rsqtmp = rsq + shift*shift + 2*sqrt(rsq)*shift;
        forceshiftfac = sqrt(rsqtmp/rsq);
        rsq = rsqtmp;
      }

      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];
      if (rsq >= params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,EFLAG,evdwl);

      // correct force for shift in rsq

      if (SHIFT_FLAG) fpair *= forceshiftfac;

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (EVFLAG) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // if mask is present, do not compute those three body terms
    if (gflag && gflag[i]) continue;

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;

    for (jj = 0; jj < numshort; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      if (SHIFT_FLAG)
        rsq1 += shift*shift + 2*sqrt(rsq1)*shift;

      if (rsq1 >= params[iparam_ij].cutsq) continue;

      const double r1inv = 1.0/sqrt(dot3(delr1, delr1));
      scale3(r1inv, delr1, r1_hat);

      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (SHIFT_FLAG)
          rsq2 += shift*shift + 2*sqrt(rsq2)*shift;

        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        double r2inv = 1.0/sqrt(dot3(delr2, delr2));
        scale3(r2inv, delr2, r2_hat);

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,r1_hat,r2_hat);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fforce,prefactor,EFLAG,evdwl);

      fpair = fforce*r1inv;

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;

      if (EVFLAG) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (SHIFT_FLAG)
          rsq2 += shift*shift + 2*sqrt(rsq2)*shift;

        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        double r2inv = 1.0/sqrt(dot3(delr2, delr2));
        scale3(r2inv, delr2, r2_hat);

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,r1_hat,r2_hat,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (VFLAG_ATOM) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

