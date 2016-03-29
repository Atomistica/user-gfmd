#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <mpi.h>

#include "li_berger.h"
#include "linearalgebra.h"

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

  // These were the default values
  Poisson_number_ = 0.5;
  shear_modulus_ = 1.0;

  if (*carg >= narg+1) {
    error_->all(FLERR,"NonperiodicStiffnessKernel::NonperiodicStiffnessKernel: "
		"Expected Poisson number and shear modulus parameters.");
  }

  Poisson_number_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"NonperiodicStiffnessKernel::NonperiodicStiffnessKernel: "
		"Can't convert nu to float.");
  }
  (*carg)++;

  shear_modulus_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"NonperiodicStiffnessKernel::NonperiodicStiffnessKernel: "
		"Can't convert myu to float.");
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
  if (complex_buffer_ && ( nx_ != nx || ny_ != ny)) {
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
        sq_uniform_PQ(shear_modulus_,Poisson_number_, 0,fac,0,
                      dx2*reglen,dy2*reglen, i*dx,j*dy,
                      MEL(3,G,0,0),MEL(3,G,1,0),MEL(3,G,2,0));
        sq_uniform_PQ(shear_modulus_,Poisson_number_, 0,0,fac,
                      dx2*reglen,dy2*reglen, i*dx,j*dy,
                      MEL(3,G,0,1),MEL(3,G,1,1),MEL(3,G,2,1));
        sq_uniform_PQ(shear_modulus_,Poisson_number_, fac,0,0,
                      dx2*reglen,dy2*reglen, i*dx,j*dy,
                      MEL(3,G,0,2),MEL(3,G,1,2),MEL(3,G,2,2));

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

#if 0 // Debug files out
    char dump_prefix[1024], fn[1024];
    for (int iiiii=0; iiiii<9; iiiii++) {
      sprintf(fn, "myGF.%i.out", iiiii);
      FILE *fuP = fopen(fn, "w");
      for (int jjjj=0; jjjj<ny; jjjj++) {
	for (int iiii=0; iiii<nx; iiii++) fprintf(fuP, "%f ", complex_buffer_[iiiii][jjjj*nx+iiii][0]);  fprintf(fuP, "\n"); 
      } //TAS // IDL  c=read_ascii("myGF.8.out") & d=c.field01 & tvscl, congrid(shift(d,16,16),128,128)
      fclose(fuP);  
    }
#endif

#ifndef FFT_FFTW3
#error NonperiodicStiffnessKernel requires FFTW3
#endif
   
    /* Fourier transform each matrix entry */
    for (int a = 0; a < 9; a++) {
      fftw_execute_dft(plan, complex_buffer_[a], complex_buffer_[a]);
    }
    
    /* Destroy the FFT plan (and release buffers) */
    fftw_destroy_plan(plan);

#if 0 /* over write the moderate and high q components with periodic GF result */ // TAS  
    
    // This section may need to be upgraded to allow negative 'fac' // TAS
    double myu_ = shear_modulus_;
    double alpha = 1.0/(2.0*(1.0-Poisson_number_));
    for (int ii = 0; ii < nx; ii++) {
      for (int jj = 0; jj < ny; jj++) {

        int i = ii, j = jj;
        i = (i <= nx2) ? i : i-nx;
        j = (j <= ny2) ? j : j-ny;

	double qx = 2.0*M_PI*i/nx;
	double qy = 2.0*M_PI*j/ny;
        double beta = sqrt(pow(qx,2.0)+pow(qy,2.0));

	int q_index = ii*ny+jj;

        if ((abs(i) > 2) || (abs(j) > 2)) {

	  complex_buffer_[0*3+0][q_index][0] = 
          1.0/(myu_*beta) - 
	    (pow(qx,2.0)*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0)));

	  complex_buffer_[1*3+1][q_index][0] = 
	    (1.0/(myu_*beta) -
	     (pow(qy,2.0)*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0))));

	  complex_buffer_[2*3+2][q_index][0] = 
	     (1.0/(2.0*myu_*alpha*beta));

	  complex_buffer_[0*3+1][q_index][0] = 
	     -qx*qy*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0));

	  complex_buffer_[1*3+0][q_index][0] = 
	     -qx*qy*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0));

	  complex_buffer_[0*3+2][q_index][1] = 
	     -qx*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0));

	  complex_buffer_[1*3+2][q_index][1] = 
	     -qy*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0));

	  complex_buffer_[2*3+0][q_index][1] = 
	     qx*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0));

	  complex_buffer_[2*3+1][q_index][1] = 
 	     qy*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0));
	}
      }
    }
#endif

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
  error_->all(FLERR,"get_stiffness_matrix without grid size not supported.");
}


void NonperiodicStiffnessKernel::get_stiffness_matrix(int nx, int ny,
                                                      double qx, double qy,
                                                      double_complex *phi,
                                                      double_complex dU)
{
  if (dU != 0.0) {
    error_->all(FLERR,"dU != 0 not supported.");
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
