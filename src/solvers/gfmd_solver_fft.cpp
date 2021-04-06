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
#include <string.h>
#include <stdlib.h>

#include "gfmd_solver_fft.h"

#include "pointers.h"

#include "comm.h"
#include "domain.h"
#include "linearalgebra.h"
#include "memory.h"
#include "mpi.h"
#include "fft3d_wrap.h"

#include "gfmd_misc.h"

using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

/* ----------------------------------------------------------------------
 * Solver superclass using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

GFMDSolverFFT::GFMDSolverFFT(LAMMPS *lmp) : GFMDSolver(lmp)
{
  strcpy(name, "GFMDSolverFFT");

  fft = NULL;
  fft_data = NULL;
}


GFMDSolverFFT::~GFMDSolverFFT()
{
  if (fft_data)
    memory->sfree(fft_data);

  if (fft)
    delete fft;
}


void GFMDSolverFFT::set_grid_size(int in_nx, int in_ny, int in_dof)
{
  GFMDSolver::set_grid_size(in_nx, in_ny, in_dof);

  int nfft;

  fft  = new FFT3d(lmp, world,
               //  fast, med,              slow
                   1,    ny,               nx,
               //
                   0, 0, ylo_loc, yhi_loc, xlo_loc, xhi_loc,
                   0, 0, ylo_loc, yhi_loc, xlo_loc, xhi_loc,
                   0, 0, &nfft,
               //  usecollective which refers to something with MPI
                   0);

  /*
  if (nfft != nxy_loc) {
    char errstr[1024];
    sprintf(errstr, "Local buffer size (= %zu) does not match buffer size "
            "return by the FFT wrapper (= %i).", nxy_loc, nfft);
    error->one(FLERR,errstr);
  }
  */

  fft_data = (double *) memory->smalloc(nxy_loc*2*sizeof(double),
                                        "GFMDSolver::fft_data");
}


void GFMDSolverFFT::fft_forward(double **input_buffer,
                                double_complex **q_buffer,
                                double_complex **fac)
{
  for (int idim = 0; idim < ndof; idim++) {
    /*
     * Fill FFT buffer and perform FFT on local data
     */
    int m = 0;
    for (int idx = 0; idx < nxy_loc; idx++) {
      /*
       * displacement U is all real.  space it out into everyother element
       */
      fft_data[m++] = input_buffer[idim][idx];
      fft_data[m++] = 0.;
    }

    /*
     * perform the FFT!
     */
    fft->compute(fft_data, fft_data, 1);

    if (fac) {
      /*
       * Multiply with phase factor
       */
      double_complex *cfac = fac[idim/3];
      m = 0;
      for (int idq = 0; idq < nxy_loc; idq++){
        /*
         * pack double array data into complex array
         */
        q_buffer[idq][idim] = 
          cfac[idq]*COMPLEX_NUMBER(fft_data[m], fft_data[m+1]);
        m += 2;
      }
    }
    else {
      /*
       * No multiplication with a phase factor
       */
      m = 0;
      for (int idq = 0; idq < nxy_loc; idq++){
        /*
         * pack double array data into complex array
         */
        q_buffer[idq][idim] = COMPLEX_NUMBER(fft_data[m], fft_data[m+1]);
        m += 2;
      }
    }
  }
}


void GFMDSolverFFT::fft_reverse(double_complex **q_buffer,
                                double **output_buffer,
                                double_complex **fac)
{
  for (int idim = 0; idim < ndof; idim++){
    if (fac) {
      /*
       * Now transform F(q) to F(r), multiply with phase factor and fill buffer
       */
      double_complex *cfac = fac[idim/3];
      int m = 0;
      for (int idq = 0; idq < nxy_loc; idq++){
        double_complex v = conj(cfac[idq])*q_buffer[idq][idim];
        fft_data[m++] = creal(v);
        fft_data[m++] = cimag(v);
      }
    }
    else {
      /*
       * No multiplication with a phase factor
       */
      int m = 0;
      for (int idq = 0; idq < nxy_loc; idq++){
        fft_data[m++] = creal(q_buffer[idq][idim]);
        fft_data[m++] = cimag(q_buffer[idq][idim]);
      }
    }

    /*
     * perform the FFT!
     */
    fft->compute(fft_data, fft_data, -1);

    /*
     * Copy buffer to output
     */
    int m = 0;
    for (int idx = 0; idx < nxy_loc; idx++) {
      /*
       * displacement U is all real.  space it out into everyother element
       */
      output_buffer[idim][idx] = fft_data[m++];
      m++;
    }
  }
}


double GFMDSolverFFT::memory_usage()
{
  double bytes = 0.0;

  // fft_data
  bytes += nxy_loc*2*sizeof(double);

  return bytes;
}


void GFMDSolverFFT::dump(char *dump_prefix, double_complex **q_buffer)
{
  if (nprocs > 1)
    error->all(FLERR,"Can only dump from single processor run.");

  char fn[1024];
#ifdef HAVE_C99
  FILE *fur[ndof], *ffr[ndof], *fui[ndof], *ffi[ndof];
  double_complex F_q[ndof];
#else
  FILE *fur[MAX_NDOF], *ffr[MAX_NDOF], *fui[MAX_NDOF], *ffi[MAX_NDOF];
  double_complex *F_q = new double_complex[ndof];
#endif

  for (int idof = 0; idof < ndof; idof++) {
    sprintf(fn, "%s.q.u%i.real.out", dump_prefix, idof);
    fur[idof] = fopen(fn, "w");
    sprintf(fn, "%s.q.u%i.imag.out", dump_prefix, idof);
    fui[idof] = fopen(fn, "w");
    sprintf(fn, "%s.q.f%i.real.out", dump_prefix, idof);
    ffr[idof] = fopen(fn, "w");
    sprintf(fn, "%s.q.f%i.imag.out", dump_prefix, idof);
    ffi[idof] = fopen(fn, "w");
  }
  sprintf(fn, "%s.q.uP.out", dump_prefix);
  FILE *fuP = fopen(fn, "w");
  sprintf(fn, "%s.q.fP.out", dump_prefix);
  FILE *ffP = fopen(fn, "w");
  sprintf(fn, "%s.q.e.out", dump_prefix);
  FILE *fe = fopen(fn, "w");

  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      int iloc = ix*ny_loc + iy;

      MatMulVec(ndof, phi[iloc], q_buffer[iloc], F_q);

      double_complex uP = 0.0, fP = 0.0, e = 0.0;
      for (int idof = 0; idof < ndof; idof++) {
        fprintf(fur[idof], " %20.10e ", creal(q_buffer[iloc][idof]));
        fprintf(fui[idof], " %20.10e ", cimag(q_buffer[iloc][idof]));

        fprintf(ffr[idof], " %20.10e ", creal(F_q[idof]));
        fprintf(ffi[idof], " %20.10e ", cimag(F_q[idof]));

        uP += q_buffer[iloc][idof]*conj(q_buffer[iloc][idof]);
        fP += F_q[idof]*conj(F_q[idof]);

        e  += q_buffer[iloc][idof]*conj(F_q[idof]);
      }
      fprintf(fuP, " %20.10e ", creal(uP));
      fprintf(ffP, " %20.10e ", creal(fP));
      fprintf(fe, " %20.10e ", creal(e));
    }
    for (int idof = 0; idof < ndof; idof++) {
      fputc('\n', fur[idof]);
      fputc('\n', fui[idof]);
      fputc('\n', ffr[idof]);
      fputc('\n', ffi[idof]);
    }
    fputc('\n', fuP);
    fputc('\n', ffP);
    fputc('\n', fe);
  }

  for (int idof = 0; idof < ndof; idof++) {
    fclose(fur[idof]);
    fclose(fui[idof]);
    fclose(ffr[idof]);
    fclose(ffi[idof]);
  }
  fclose(fuP);
  fclose(ffP);
  fclose(fe);

#ifndef HAVE_C99
  delete [] F_q;
#endif
}


/* ----------------------------------------------------------------------
 * dump stiffness coefficients for plotting in gnuplot
 * --------------------------------------------------------------------*/

void GFMDSolverFFT::dump_stiffness()
{
  if (nprocs > 1)
    error->all(FLERR,"fix gfmd/static: Dump stiffness only works from a single "
               "processor run.");

  if (me == 0) {
#ifdef HAVE_C99
    FILE *freal[ndof][ndof], *fimag[ndof][ndof];
#else
    FILE *freal[MAX_NDOF][MAX_NDOF], *fimag[MAX_NDOF][MAX_NDOF];
#endif
    FILE *ftrreal, *ftrimag;

    for (int idim = 0; idim < ndof; idim++) {
      for (int jdim = 0; jdim < ndof; jdim++) {
        char fn[1024];
        sprintf(fn, "phi%i%i.real.out", idim, jdim);
        freal[idim][jdim] = fopen(fn, "w");
        sprintf(fn, "phi%i%i.imag.out", idim, jdim);
        fimag[idim][jdim] = fopen(fn, "w");
      }
    }
    ftrreal = fopen("phitr.real.out", "w");
    ftrimag = fopen("phitr.imag.out", "w");

    int n = 0;
    for (int j = 0; j < ny; j++) {
      for (int i = xlo_loc; i <= xhi_loc; i++) {
        int ndim = 0;
        double_complex trace = 0.0;
        for (int idim = 0; idim < ndof; idim++) {
          for (int jdim = 0; jdim < ndof; jdim++) {
            fprintf(freal[idim][jdim], " %e ", creal(phi[n][ndim]));
            fprintf(fimag[idim][jdim], " %e ", cimag(phi[n][ndim]));
            trace += phi[n][ndim];
            ndim++;
          }
        }
        n++;
      }

      for (int idim = 0; idim < ndof; idim++) {
        for (int jdim = 0; jdim < ndof; jdim++) {
          fputc('\n', freal[idim][jdim]);
          fputc('\n', fimag[idim][jdim]);
        }
      }
      fputc('\n', ftrreal);
      fputc('\n', ftrimag);
    }

    for (int idim = 0; idim < ndof; idim++) {
      for (int jdim = 0; jdim < ndof; jdim++) {
        fclose(freal[idim][jdim]);
        fclose(fimag[idim][jdim]);
      }
    }
    fclose(ftrreal);
    fclose(ftrimag);
  }
}


/* ----------------------------------------------------------------------
 * dump stiffness coefficients for plotting in gnuplot
 * --------------------------------------------------------------------*/

void GFMDSolverFFT::dump_greens_function()
{
  if (nprocs > 1)
    error->all(FLERR,"fix gfmd/static: Dump Green's function only works from a "
               "single processor run.");

  if (me == 0) {
#ifdef HAVE_C99
    FILE *freal[ndof][ndof], *fimag[ndof][ndof];
    double_complex G[ndof*ndof];
#else
    FILE *freal[MAX_NDOF][MAX_NDOF], *fimag[MAX_NDOF][MAX_NDOF];
    double_complex *G = new double_complex[ndof*ndof];
#endif
    FILE *ftrreal, *ftrimag;

    for (int idim = 0; idim < ndof; idim++) {
      for (int jdim = 0; jdim < ndof; jdim++) {
        char fn[1024];
        sprintf(fn, "g%i%i.real.out", idim, jdim);
        freal[idim][jdim] = fopen(fn, "w");
        sprintf(fn, "g%i%i.imag.out", idim, jdim);
        fimag[idim][jdim] = fopen(fn, "w");
      }
    }
    ftrreal = fopen("gtr.real.out", "w");
    ftrimag = fopen("gtr.imag.out", "w");

    int n = 0;
    for (int j = 0; j < ny; j++) {
      for (int i = xlo_loc; i <= xhi_loc; i++) {
        int ndim = 0;
        double_complex trace = 0.0;
        memcpy(G, phi[n], ndof_sq*sizeof(double_complex));
        GaussJordan(ndof, G, error);
        for (int idim = 0; idim < ndof; idim++) {
          for (int jdim = 0; jdim < ndof; jdim++) {
            fprintf(freal[idim][jdim], " %e ", creal(G[ndim]));
            fprintf(fimag[idim][jdim], " %e ", cimag(G[ndim]));
            trace += G[ndim];
            ndim++;
          }
        }
        n++;
      }

      for (int idim = 0; idim < ndof; idim++) {
        for (int jdim = 0; jdim < ndof; jdim++) {
          fputc('\n', freal[idim][jdim]);
          fputc('\n', fimag[idim][jdim]);
        }
      }
      fputc('\n', ftrreal);
      fputc('\n', ftrimag);
    }

    for (int idim = 0; idim < ndof; idim++) {
      for (int jdim = 0; jdim < ndof; jdim++) {
        fclose(freal[idim][jdim]);
        fclose(fimag[idim][jdim]);
      }
    }
    fclose(ftrreal);
    fclose(ftrimag);

#ifndef HAVE_C99
    delete [] G;
#endif
  }
}
