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
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <mpi.h>

#include "li_berger.h"
#include "linearalgebra.h"

#include "crystal_surface.h"

#include "nonperiodic_stiffness.h"

using namespace LiBerger;

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Nonperiodic elasticity kernel, i.e. open boundary conditions
 * --------------------------------------------------------------------*/

NonperiodicStiffnessKernel::NonperiodicStiffnessKernel(int narg, int *carg,
                                                       char **arg,
                                                       Domain *domain,
                                                       Memory *memory,
                                                       Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
    char *endptr;

    dim_ = 3;
    strcpy(name_, "nonperiodic");

    // These are the default values
    Poisson_number_ = 0.5;
    shear_modulus_ = 1.0;

#if 0
    int i = *carg;
    if (i >= narg)
        error->all(FLERR, "Missing crystal surface.");

    (*carg)++;
    crystal_surface_ = crystal_surface_factory(arg[i], narg, carg, arg, error);
    if (!crystal_surface_) {
        char errstr[1024];
        sprintf(errstr, "Could not find crystal surface '%s'.", arg[i]);
        error->all(FLERR,errstr);
    }
#endif

    if (*carg >= narg)
        error->all(FLERR, "Missing Poisson number.");

    Poisson_number_ = strtod(arg[*carg], &endptr);
    if (endptr == arg[*carg]) {
        char errstr[1024];
        sprintf(errstr, "Can't convert Poisson number '%s' to floating point.",
                arg[*carg]);
        error_->all(FLERR, errstr);
    }
    (*carg)++;

    if (*carg >= narg)
        error->all(FLERR, "Missing shear modulus.");

    shear_modulus_ = strtod(arg[*carg], &endptr);
    if (endptr == arg[*carg]) {
        char errstr[1024];
        sprintf(errstr, "Can't convert shear modulus '%s' to floating point.",
                arg[*carg]);
        error_->all(FLERR, errstr);
    }
    (*carg)++;

    // Can't really specify height in layers
    height_ = 1;

    // Initialize buffers to zero
    nx_ = ny_ = 0;
    complex_buffer_ = NULL;
}


NonperiodicStiffnessKernel::~NonperiodicStiffnessKernel()
{
    if (complex_buffer_) {
        memory_->destroy(complex_buffer_);
    }
}


void NonperiodicStiffnessKernel::create_and_fill_buffer(int nx, int ny)
{
    if (complex_buffer_ && (nx_ != nx || ny_ != ny)) {
        memory_->destroy(complex_buffer_);
    
        complex_buffer_ = NULL;
        
        nx_ = nx;
        ny_ = ny;
    }
  
    if (!complex_buffer_) {
        memory_->create(complex_buffer_, 9, nx*ny,
                        "NonperiodicStiffnessKernel::complex_buffer_");
        
        /* Plan the FFT */
        fftw_plan plan = fftw_plan_dft_2d(nx, ny, complex_buffer_[0],
                                          complex_buffer_[0], 1,
                                          FFTW_ESTIMATE);
    
        memset(complex_buffer_[0], 0, 9*nx*ny*sizeof(fftw_complex));
  
        double xprd = domain_->xprd;
        double yprd = domain_->yprd;
    
        double dx = xprd/nx, dy = yprd/ny;
        double dx2 = dx/2, dy2 = dy/2;
        double reglen = 1.0; // regularization lenth scale.
                                            // 0.9328 leads to match of periodic run, single displaced atom. indpendent of nu
        double fac = 1./(dx*dy*reglen*reglen); // units of force / area
        
        int nx2 = nx/2, ny2 = ny/2;
    
        bool nx_is_even = nx % 2 == 0, ny_is_even = ny % 2 == 0;
     
        int k = 0;
        for (int ii = 0; ii < nx; ii++) {
            for (int jj = 0; jj < ny; jj++, k++) {
                double G[9];
                memset(G, 0, 9*sizeof(double));
        
                int i = ii, j = jj;
                i = (i <= nx2) ? i : i-nx;
                j = (j <= ny2) ? j : j-ny;
               
                /* Get the Green's function. Note, sq_uniform_PQ takes P (pressure in
                   z-direction first). */
                sq_uniform_PQ(shear_modulus_, Poisson_number_,  0,fac,0,
                              dx2*reglen, dy2*reglen,  i*dx, j*dy,
                              MEL(3,G,0,0), MEL(3,G,1,0), MEL(3,G,2,0));
                sq_uniform_PQ(shear_modulus_, Poisson_number_,  0,0,fac,
                              dx2*reglen, dy2*reglen,  i*dx, j*dy,
                              MEL(3,G,0,1), MEL(3,G,1,1), MEL(3,G,2,1));
                sq_uniform_PQ(shear_modulus_, Poisson_number_,  fac,0,0,
                              dx2*reglen, dy2*reglen,  i*dx, j*dy,
                              MEL(3,G,0,2), MEL(3,G,1,2), MEL(3,G,2,2));

                // I verify that this gives the Boussinesq solution as reglen goes to zero

#if 1 // Setting real space furthest edge to zero displacement allows the stiffness to be Hermetian phi31 cc=phi13
                if ((nx_is_even && ii == nx2) || (ny_is_even && jj == ny2)) {
                    //MEL(3,G,0,0) = 0.0;
                    //MEL(3,G,1,1) = 0.0;
                    //MEL(3,G,2,2) = 0.0;
                    MEL(3,G,0,1) = 0.0;
                    MEL(3,G,0,2) = 0.0;
                    MEL(3,G,1,2) = 0.0;
                    MEL(3,G,1,0) = 0.0;
                    MEL(3,G,2,0) = 0.0;
                    MEL(3,G,2,1) = 0.0;
                }
#endif
                      
                /* Copy to Green's function buffer */
                for (int a = 0; a < 9; a++) {
                  complex_buffer_[a][k][0] = G[a];
                  complex_buffer_[a][k][1] = 0.0;
                }
          
            }
        }
  
  #ifndef FFT_FFTW3
  #error NonperiodicStiffnessKernel requires FFTW3
  #endif
   
        /* Fourier transform each matrix entry */
        for (int a = 0; a < 9; a++) {
            fftw_execute_dft(plan, complex_buffer_[a], complex_buffer_[a]);
        }
        
        /* Destroy the FFT plan (and release buffers) */
        fftw_destroy_plan(plan);
  
        /* Invert Green's function to obtain stiffness matrix coefficients */
        for (int i = 0; i < nx*ny; i++) {
            double_complex G[9];
          
            for (int a = 0; a < 9; a++) {
                G[a] = COMPLEX_NUMBER(complex_buffer_[a][i][0],
                                      complex_buffer_[a][i][1]);
            }
    
            invert3x3(G);
    
            for (int a = 0; a < 9; a++) {
                complex_buffer_[a][i][0] = creal(G[a]);
                complex_buffer_[a][i][1] = cimag(G[a]);
            }
        }
    }  
}


void NonperiodicStiffnessKernel::get_stiffness_matrix(double qx, double qy,
                                                      double_complex *phi,
                                                      double_complex dU)
{
    error_->all(FLERR, "get_stiffness_matrix without grid size not supported.");
}


void NonperiodicStiffnessKernel::get_stiffness_matrix(int nx, int ny,
                                                      double qx, double qy,
                                                      double_complex *phi,
                                                      double_complex dU)
{
    if (dU != 0.0) {
        error_->all(FLERR, " dU != 0 not supported.");
    }

    /* Check if buffer size has changed or this is first call */
    create_and_fill_buffer(nx, ny);

    /* Extract grid index from qx, qy */
    double x = qx*nx/(2.0*M_PI), y = qy*ny/(2.0*M_PI);
    
    int ix = round(x), iy = round(y);
    if (abs(x-ix) > 1e6 || abs(y-iy) > 1e6) {
        error_->one(FLERR,"Could not extract grid index from q-vector.");
    }
    
    /* Wrap into range */
    if (ix < 0)  ix += nx;
    if (iy < 0)  iy += ny;
   
    /* Copy from buffer */
    int i = ix*ny+iy;
    for (int a = 0; a < 9; a++) {
        phi[a] = COMPLEX_NUMBER(complex_buffer_[a][i][0],
                                complex_buffer_[a][i][1]);
    }
    
    //printf("ix = %i, iy = %i\n", ix, iy);
    //printmat(3, phi);
}


} /* namespace */
