/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

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
#include "gfmd_misc.h"

#include "linearalgebra.h"

/* ----------------------------------------------------------------------
 * Loop over all wavevectors qx, qy and call the stiffness kernel to
 * get the appropriate stiffness matrix.
 * --------------------------------------------------------------------*/

void fill_phi_buffer(int ndof,
                     int nx, int xlo_loc, int xhi_loc,
                     int ny, int ylo_loc, int yhi_loc,
                     StiffnessKernel *kernel, double_complex **phi_q,
                     bool normalize, Error *error)
{ // printf("Called fill_phi_buffer\n"); // TAS
    kernel->pre_compute();

    int nu = ndof/3;
    int idq = 0;
    for (int i = xlo_loc; i <= xhi_loc; i++) {
        double qx = (i <= int((nx)/2)) ?
            (2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
        for (int j = ylo_loc; j <= yhi_loc; j++) {
            double qy = (j <= int((ny)/2)) ? 
                (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

            // This is where our GFunctions^-1 are imported from
            // surface_stiffness.
            kernel->get_stiffness_matrix(nx, ny, qx, qy, phi_q[idq]);

            // Sanity check 1: Check if matrix is Hermitian.
#if 1
            // If compliance is ever infinite in a certain direction, this will
            // false positive.
            if (!_is_Hermitian(ndof, phi_q[idq], 1e-6)) { // tol for debugging
                char errstr[1024];
                sprintf(errstr, "Kernel is not Hermitian at q = 2 pi "
                        "(%f, %f)!\n", qx/(2.0*M_PI), qy/(2.0*M_PI));

                printf("here is the un hermitian matrix:\n");
                printmat(ndof,phi_q[idq]);
                printf(" nu %d ndof %d\n",nu, ndof);
                double_complex myout[144];
                conj_transpose(ndof, myout, phi_q[idq]);
                printf("here is the cctrans: \n");
                printmat(ndof, myout);
                for (int i=0; i<ndof*ndof; i++) myout[i] -= (phi_q[idq])[i];
                printf("here is the difference which should be 0\n");
                printmat(ndof, myout);
                error->all(FLERR,errstr);
            }
#endif

            // Sanity check 2: Check if matrix is positive definite.
            // (FIXME!!! Requires LAPACK.)
#if GFMD_CHECK_POSITIVE_DEFINITE
            char name[1024];
            sprintf(name, "phi(%f,%f)", qx/(2.0*M_PI), qy/(2.0*M_PI));
            warn_positive_definite(name, ndof, phi_q[idq], error);
#endif

            idq++;
        }
    }

    if (normalize) {
        int nxy_loc = (xhi_loc-xlo_loc+1)*(yhi_loc-ylo_loc+1);
        int ndof_sq = ndof*ndof;
        // Divide phi_q by (nx*ny) to reduce float operation after FFT.
        double inv_nxny = 1.0/double(nx*ny);
        for (int idq = 0; idq < nxy_loc; idq++) {
            for (int idim = 0; idim < ndof_sq; idim++) {
                // Premptively divide for later conven.
                phi_q[idq][idim] *= inv_nxny;
            }
        }
    }

    kernel->post_compute();
}