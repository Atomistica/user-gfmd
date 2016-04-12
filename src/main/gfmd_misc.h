/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */
#ifndef __GFMD_MISC_H
#define __GFMD_MISC_H

#include "surface_stiffness.h"


void fill_phi_buffer(int, int, int, int, int, int, int, StiffnessKernel *,
                     double_complex **, bool, Error *);

/* ----------------------------------------------------------------------
 * Precondition gradient g
 * --------------------------------------------------------------------*/

#define DEF_G 1
#define TR_G 2

template<int t>
void precondition_gradient(int ndof, int nq, double *cavg,
                           double_complex **phi,
                           double_complex **g,
                           Error *error)
{
  /*
   * Perform matrix operation: gP(q) = (Phi(q) + Cavg)^(-1) x g(q)
   * and compute potential energy
   */

  if (ndof == 1) {

    // Allow for loop unrolling in 1 dimensions

#pragma omp parallel for default(none) shared(phi, g) firstprivate(cavg, nq)
    for (int idq = 0; idq < nq; idq++) {
    double_complex M;
      M = 1.0/(phi[idq][0] + cavg[0]);
      if (t == DEF_G) {
        g[idq][0] = M*g[idq][0];
      }
      else {
        g[0][idq] = M*g[0][idq];
      }
    }

  }
  else if (ndof == 3) {

    // Allow for loop unrolling in 3 dimensions

#pragma omp parallel for default(none) shared(phi, g) firstprivate(cavg, nq)
    for (int idq = 0; idq < nq; idq++) {
      double_complex M[9], F_q[3];
      for (int idim = 0; idim < 9; idim++) {
        M[idim] = phi[idq][idim];
      }
      for (int i = 0; i < 9; i++) {
        M[i] += cavg[i];
      }
      invert3x3(M);
      if (t == DEF_G) {
        mat_mul_vec(3, M, g[idq], F_q);
        for (int idim = 0; idim < 3; idim++) {
          g[idq][idim] = F_q[idim];
        }
      }
      else {
        double_complex gi[3];
        for (int idim = 0; idim < 3; idim++) {
          gi[idim] = g[idim][idq];
        }
        mat_mul_vec(3, M, gi, F_q);
        for (int idim = 0; idim < 3; idim++) {
          g[idim][idq] = F_q[idim];
        }
      }
    }

  } 
  else {

    int ndof_sq = ndof*ndof;
#pragma omp parallel for default(none) shared(error, phi, g)    \
  firstprivate(cavg, ndof, ndof_sq, nq)
    for (int idq = 0; idq < nq; idq++) {
      double_complex M[ndof_sq], F_q[ndof];
      for (int idim = 0; idim < ndof_sq; idim++) {
        M[idim] = phi[idq][idim];
      }
      for (int i = 0; i < ndof_sq; i++) {
        M[i] += cavg[i];
      }
      GaussJordan(ndof, M, error);
      if (t == DEF_G) {
        mat_mul_vec(ndof, M, g[idq], F_q);
        for (int idim = 0; idim < 3; idim++) {
          g[idq][idim] = F_q[idim];
        }
      }
      else {
        double_complex gi[ndof];
        for (int idim = 0; idim < ndof; idim++) {
          gi[idim] = g[idim][idq];
        }
        mat_mul_vec(ndof, M, gi, F_q);
        for (int idim = 0; idim < ndof; idim++) {
          g[idim][idq] = F_q[idim];
        }
      }
    }

  }
}


#endif
