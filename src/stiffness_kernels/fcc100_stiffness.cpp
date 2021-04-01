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
// Eqn labels refer to Pastewka et. al. PRB 86 075459 (2012)

// Ex: For third-nearest neighbors, LJ,s=1,e=1
//IDL> a= 1.10069d
//IDL> b= 1.5566056d
//IDL> c= 1.90644d
//IDL> print, -48d / a^(13d)  +24d /a^(7d)
//      -1.5291539
//IDL> print, -48d / b^(13d)  +24d /b^(7d)
//      0.93145283
//IDL> print, -48d / c^(13d)  +24d /c^(7d)
//      0.25128704


#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "fcc100_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Face centered cubic (100) surface
 * --------------------------------------------------------------------*/

FCC100StiffnessKernel::FCC100StiffnessKernel(int narg, int *carg, char **arg,
					     Domain *domain, Memory *memory,
					     Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  strcpy(name_, "fcc100");

  z_ = 1;
  mulphase_ = 0.5;

  // Input order: k1 (kappa1) k2 (kappa2) k3 (kappa3)
  k_[0] = 0.0;
  k_[1] = 0.0;
  k_[2] = 0.0;
  kappa_[0] = 0.0;
  kappa_[1] = 0.0;
  kappa_[2] = 0.0;

  if (*carg >= narg) {
    error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Expected "
		"number of interacting neighbors.");
  }

  //int l = strlen(arg[*carg]);
  //bool x = false;
  //// x denotes kappa to be set. 
  //if (arg[*carg][l-1] == 'x') { // TAS2
  //  printf("Setting the x option for kappa\n");
  //  x = true;
  //  arg[*carg][l] = '\0';
  //}

  // Before the the nu_ input location, find optional key word 'x' (ending in x)
  // Choosing to set kappa values
  bool x = false;
  int l = strlen(arg[*carg]);
  if (arg[*carg][l-1] == 'x') {
    x = true;
    printf("Setting the x option for kappa\n");
    (*carg)++;
  }

  // Number of layers: nu_
  nu_ = strtol(arg[*carg], &endptr, 10); 
  printf("nu_ is integer %s \n",arg[*carg]);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Error "
		"converting number of layers to integer.");
  }
  (*carg)++;

  if (nu_ < 1)
    error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Number "
		"of bounday GF layers must be at least 1.");

  // Number of spring constants defined for the interacting neighbors: nk_
  switch( nu_ )
  {
    case 1: nk_ = 1; break; // 1 layer means only support up to nn in the fcc lattice
    case 2: nk_ = 3; break; // 2 layers can support nn, 2n, and 3n spring constants
    default: nk_ = 3; break;
  }

  printf("Expecting %d spring constants.  They can be zero.\n", nk_);
  //char errstr[1024];
  //sprintf(errstr, "Expecting %i spring constants.  They can be zero.", nk_);
  //error_->all(FLERR,errstr);

  for (int i = 0; i < nk_; i++) {
    if (*carg >= narg) {
      error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Spring "
		  "constant expected.");
    }

    printf("k_ is the float: %s \n",arg[*carg]);
    k_[i] = strtod(arg[*carg], &endptr);
    if (endptr == arg[*carg]) {
      error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Error "
		  "converting spring constant to number.");
    }
    (*carg)++;

    if (x) {
      if (*carg >= narg) {
	error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: "
		    "Spring constant expected.");
      }

    printf("kappa_ is the float: %s \n",arg[*carg]);
      kappa_[i] = strtod(arg[*carg], &endptr);
      if (endptr == arg[*carg]) {
	error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Error "
		    "converting spring constant to number.");
      }
      (*carg)++;
    }
  }

  linf_ = new double[nu_];
  memset(linf_, 0, nu_*sizeof(double));

  bool keyword_found = true;
  while (*carg < narg-1 && keyword_found) {
    keyword_found = false;
    if (!strcmp(arg[*carg], "height")) {
      (*carg)++;
      printf("height is the integer: %s\n",arg[*carg]);
      height_ = strtol(arg[*carg], &endptr, 10);
      if (endptr == arg[*carg]) {
	error_->all(FLERR,"FCC100StiffnessKernel::FCC100StiffnessKernel: Error "
		    "converting height argument to number.");
      }
      (*carg)++;
      keyword_found = true;
    }
    // add a linear force (per atom) contribution normal to the surface to
    // each layer to offset surface relaxation
    else if (!strcmp(arg[*carg], "linf")) {
      double linfsum = 0.0;
      for (int i = 0; i < nu_; i++) {
        (*carg)++;
        if ((*carg) >= narg) {
          char errstr[1024];
          sprintf(errstr, "fix gfmd: Stiffness kernel report %i layers,"
                  " please provide that many *linf* parameters.", nu_);
          error_->all(FLERR, errstr);
        }
        linf_[i] = atof(arg[*carg]);
        linfsum += linf_[i];
      }
      if (abs(linfsum) > 1e-12) {
	char errstr[1024];
        sprintf(errstr, "fix gfmd: Linear force contributions must sum to "
		"zero but sums to %e.", linfsum);
	error_->all(FLERR, errstr);
      }
      (*carg)++;
      keyword_found = true;
    }
    else if (!strcmp(arg[*carg], "normal")) {
      (*carg)++;
      z_ = strtol(arg[*carg], &endptr, 10);
      if (endptr == arg[*carg]) {
	char errstr[1024];
	sprintf(errstr, "fix gfmd: Could not convert surface normal (%s) to "
		"integer.", arg[*carg]);
	error_->all(FLERR, errstr);
      }

      if (z_ != 1 && z_ != -1) {
	char errstr[1024];
	sprintf(errstr, "fix gfmd: Surface normal must be either +1 or -1, but "
		"is %i.", z_);
	error_->all(FLERR, errstr);
      }

      (*carg)++;
      keyword_found = true;
    }
    else if (!strcmp(arg[*carg], "mulphase")) {
      (*carg)++;
      mulphase_ = strtod(arg[*carg], &endptr);
      printf("mulphase  %s \n",arg[*carg]);
      if (endptr == arg[*carg]) {
        char errstr[1024];
        sprintf(errstr, "fix gfmd: Could not convert mulphase (%s) to "
                "double.", arg[*carg]);
        error_->all(FLERR, errstr);
      }
      (*carg)++;
      keyword_found = true;
    }
  }
  /*
  if (logfile) {
    fprintf(logfile, "FCC100 boundary with %i layers and normal in direction "
	    "%i.", nu_, z_);
  }
  */

  // natomsunitcell_ = nu_; // For now assume one atom per layer in the unit cell
  dim_ = ndof_*nu_; // Tot dim=3*nc (in multiatom cell) or 3*num (1-atom cell) // TAS
  npars_ = 2*nk_; // TAS number of parameters, obsolete
}

FCC100StiffnessKernel::~FCC100StiffnessKernel()
{
  delete [] linf_;
}

void FCC100StiffnessKernel::set_parameters(double *pars)
{
  for (int i = 0; i < nk_; i++) {
    k_[i] = pars[i];
    kappa_[i] = pars[nk_+i];
  }
}

void FCC100StiffnessKernel::dump_info(FILE *f)
{
  StiffnessKernel::dump_info(f);

  fprintf(f, "Surface normal is %i.\n", z_);
  fprintf(f, "Phase multiplier is %f.\n", mulphase_);

  for (int i = 0; i < nk_; i++) {
    fprintf(f, "Spring constants for %i-th neighbor interaction are %f (k) and "
	    "%f (kappa).\n", i+1, k_[i], kappa_[i]);
  }

  fprintf(f, "Linear forces on each atom = ");
  for (int i = 0; i < nu_; i++)  fprintf(f, "%f ", linf_[i]);
  fprintf(f, "\n");
}

// Multi-atom unit cells written as single atom cells, per Appendix A.3
// Supports up to nu_=4 layers and interaction range 2 layers
void FCC100StiffnessKernel::get_per_layer_dynamical_matrices(double qx,
							     double qy,
							     double_complex **D,
							     double *xcshift, double *ycshift)
{
  double cx, cy, cx2, cy2, sx, sy, sx2, sy2, sqrt2;
  double_complex *U0, *U, *V, *U1, *W;
  // nu_ up to 4: U0 U1 U2 U3 U      
  // interaction up to 2: U,V,W
  // pair potential so: the following are required by D, but are slaved to the above
  double_complex *U2, *U3, *V0, *V1, *V2, *V3, *W0, *W1, *W2, *W3; 

  cx = cos(qx);
  cy = cos(qy);
  cx2 = cos(qx/2);
  cy2 = cos(qy/2);
  sx = sin(qx);
  sy = sin(qy);
  sx2 = sin(qx/2);
  sy2 = sin(qy/2);
  sqrt2 = sqrt((double)(2.0));

  // D: list of pointers to the diagonal (ie non-zero) dynam matrix elements
  // If there are nu GF boundary layers, consider them as nu single-atom cells
  // To match paper notation, nu-1 tildas sit above each variable name

  U0 = D[0];
  U  = D[nu_];
  V0 = D[nu_+1];
  V  = D[2*nu_+1];
  if (nu_ > 1) {
    U1 = D_[1];
    V1 = D_[nu_+2];
    W0 = D_[2*nu_+2];
    W1 = D_[2*nu_+3];
    W  = D_[3*nu_+2];
  }
  if (nu_ > 2) {
    U2 = D_[2];
    V2 = D_[nu_+3];
    W2 = D_[2*nu_+4];
  }
  if (nu_ > 3) {
    U3 = D_[3];
    V3 = D_[nu_+4];
    W3 = D_[2*nu_+5];
  }


  memset(D[0], 0, (nu_+1)*(nu_+1)*9*sizeof(double_complex));

  /*
   * Nearest-neighbor interactions
   */
  // arithmetic order chosen so integer/float operations promoted to double
  /* U is diagonal */ // (Eqn A12, A13)
  MEL(ndof_, U, 0, 0) = k_[0]*(4 - 2*cx) + kappa_[0]*(8 - 2*cy);
  MEL(ndof_, U, 1, 1) = k_[0]*(4 - 2*cy) + kappa_[0]*(8 - 2*cx);
  MEL(ndof_, U, 2, 2) = k_[0]*4          + kappa_[0]*(8 - 2*cx - 2*cy);

  /* U0 is diagonal */ // (Eqn A14, A15)
  MEL(ndof_, U0, 0, 0) = MEL(ndof_, U, 0, 0) - k_[0]*1 - kappa_[0]*3;
  MEL(ndof_, U0, 1, 1) = MEL(ndof_, U, 1, 1) - k_[0]*1 - kappa_[0]*3;
  MEL(ndof_, U0, 2, 2) = MEL(ndof_, U, 2, 2) - k_[0]*2 - kappa_[0]*2;

  /* diagonal of V */ // (Eqn A16, A17)
  MEL(ndof_, V, 0, 0) = -k_[0]*cx2*cy2   - kappa_[0]*3*cx2*cy2;
  MEL(ndof_, V, 1, 1) = -k_[0]*cx2*cy2   - kappa_[0]*3*cx2*cy2;
  MEL(ndof_, V, 2, 2) = -k_[0]*2*cx2*cy2 - kappa_[0]*2*cx2*cy2;

  /* off-diagonal of V */ // (Eqn A16, A17)
  MEL(ndof_, V, 0, 1) =      k_[0]*sx2*sy2       - kappa_[0]*sx2*sy2;
  MEL(ndof_, V, 0, 2) =  z_*(k_[0]*sqrt2*sx2*cy2 - kappa_[0]*sqrt2*sx2*cy2)*cI;
  MEL(ndof_, V, 1, 2) =  z_*(k_[0]*sqrt2*cx2*sy2 - kappa_[0]*sqrt2*cx2*sy2)*cI;
  
  MEL(ndof_, V, 1, 0) =      k_[0]*sx2*sy2       - kappa_[0]*sx2*sy2;
  MEL(ndof_, V, 2, 0) =  z_*(k_[0]*sqrt2*sx2*cy2 - kappa_[0]*sqrt2*sx2*cy2)*cI;
  MEL(ndof_, V, 2, 1) =  z_*(k_[0]*sqrt2*cx2*sy2 - kappa_[0]*sqrt2*cx2*sy2)*cI;

  if (nu_ >= 2) {
    /*
     * Interactions that extend vertically a full fcc lattice cell:
     * 2n (second nearest neighbor) and 3n interactions now included.
     */

    /* Now U,V,W variables have implied tildas to match PRB paper.
     * Third spring constants (k_[2]) not shown in the PRB paper.
     */ 

    double_complex U_2nd[9];

    memset(U_2nd, 0, 9*sizeof(double_complex));

    /* diagonal of U^(2) */ // (Eqn A21, A22)
#if 1
    MEL(ndof_, U_2nd, 0, 0) = k_[1]*(2-2*cx*cy) + kappa_[1]*(4-2*cx*cy)
                            + k_[2]*8           + kappa_[2]*16;
    MEL(ndof_, U_2nd, 1, 1) = k_[1]*(2-2*cx*cy) + kappa_[1]*(4-2*cx*cy)
                            + k_[2]*8           + kappa_[2]*16;
    MEL(ndof_, U_2nd, 2, 2) = k_[1]*2           + kappa_[1]*(4-4*cx*cy)
                            + k_[2]*8           + kappa_[2]*16;
#endif

    /* off-diagonal of U^(2) */ // (Eqn A21, A22)
#if 1
    MEL(ndof_, U_2nd, 0, 1) = k_[1]*2*sx*sy - kappa_[1]*2*sx*sy;
    MEL(ndof_, U_2nd, 1, 0) = k_[1]*2*sx*sy - kappa_[1]*2*sx*sy;
#endif

#if 1 
    // Utilda_0 = (Eqn A23, A25 1st term + 2nd term now add 3rd term)
    // Utilda = Utildann + Utilda2n (Eqn unnumbered after A20)
    for (int i = 0; i < ndof_; i++) {
      for (int j = 0; j < ndof_; j++) {
	MEL(ndof_, U0, i, j) += MEL(ndof_, U_2nd, i, j);
	MEL(ndof_, U,  i, j) += MEL(ndof_, U_2nd, i, j);
      }
    }
#endif

    /* U1 = U + (U^(2) - ...) */
    memcpy(U1, U, 9*sizeof(double_complex));

#if 1
    // (Eqn A23 fourth term, A25 fourth term)
    MEL(ndof_, U0, 0, 0) -= kappa_[1] + k_[2]*4 + kappa_[2]*8;
    MEL(ndof_, U0, 1, 1) -= kappa_[1] + k_[2]*4 + kappa_[2]*8;
    MEL(ndof_, U0, 2, 2) -= k_[1]     + k_[2]*4 + kappa_[2]*8;

    // (Eqn A24 fourth term, A26 fourth term)
    MEL(ndof_, U1, 0, 0) -= kappa_[1] + k_[2]*2/3 + kappa_[2]*10/3;
    MEL(ndof_, U1, 1, 1) -= kappa_[1] + k_[2]*2/3 + kappa_[2]*10/3;
    MEL(ndof_, U1, 2, 2) -= k_[1]     + k_[2]*8/3 + kappa_[2]*4/3;


    // diagonal of V
    MEL(ndof_, V, 0, 0) +=    -k_[2]*1/3*cx2*cy2*(-10+18*cx+2*cy)	
                              -kappa_[2]*1/3*cx2*cy2*(-14+6*cx+22*cy); 
    MEL(ndof_, V, 1, 1) +=    -k_[2]*1/3*cx2*cy2*(-10+2*cx+18*cy)	
                              -kappa_[2]*1/3*cx2*cy2*(-14+22*cx+6*cy);
    MEL(ndof_, V, 2, 2) +=    -k_[2]*2/3*cx2*cy2*(-2+2*(cx+cy))
                              -kappa_[2]*2/3*cx2*cy2*(-10+10*(cx+cy));
   
    // off-diagonal of V
    MEL(ndof_, V, 0, 1) +=    -k_[2]*(-2-2*(cx+cy))*sx2*sy2	
                              -kappa_[2]*(2+2*(cx+cy))*sx2*sy2;
    MEL(ndof_, V, 0, 2) +=z_*(-k_[2]*1/3*sqrt2*cy2*(-2-2*(3*cx+cy))*sx2
                              -kappa_[2]*1/3*sqrt2*cy2*(2+2*(3*cx+cy))*sx2)*cI;
    MEL(ndof_, V, 1, 2) +=z_*(-k_[2]*1/3*sqrt2*cx2*(-2-2*(cx+3*cy))*sy2
                              -kappa_[2]*1/3*sqrt2*cx2*(2+2*(cx+3*cy))*sy2)*cI;
   
    MEL(ndof_, V, 1, 0) +=    -k_[2]*(-2-2*(cx+cy))*sx2*sy2	
                              -kappa_[2]*(2+2*(cx+cy))*sx2*sy2;
    MEL(ndof_, V, 2, 0) +=z_*(-k_[2]*1/3*sqrt2*cy2*(-2-2*(3*cx+cy))*sx2
                              -kappa_[2]*1/3*sqrt2*cy2*(2+2*(3*cx+cy))*sx2)*cI;
    MEL(ndof_, V, 2, 1) +=z_*(-k_[2]*1/3*sqrt2*cx2*(-2-2*(cx+3*cy))*sy2
                              -kappa_[2]*1/3*sqrt2*cx2*(2+2*(cx+3*cy))*sy2)*cI;

    // (Eqn A29, A30)
    MEL(ndof_, W, 0, 0) = -kappa_[1] 
                         - k_[2]*2/3*cx      -kappa_[2]*(4*cx/3 + 2*cy);
    MEL(ndof_, W, 1, 1) = -kappa_[1] 
                         - k_[2]*2/3*cy      -kappa_[2]*(2*cx + 4*cy/3);
    MEL(ndof_, W, 2, 2) = -k_[1]     
                         - k_[2]*4/3*(cx+cy) -kappa_[2]*2/3*(cx + cy);
    // off-diagonal of W
    MEL(ndof_, W, 0, 1) = 0;
    MEL(ndof_, W, 0, 2) = z_*(k_[2]*2/3*sqrt2*sx - kappa_[2]*2/3*sqrt2*sx)*cI;
    MEL(ndof_, W, 1, 2) = z_*(k_[2]*2/3*sqrt2*sy - kappa_[2]*2/3*sqrt2*sy)*cI;

    MEL(ndof_, W, 1, 0) = 0;
    MEL(ndof_, W, 2, 0) = z_*(k_[2]*2/3*sqrt2*sx - kappa_[2]*2/3*sqrt2*sx)*cI;
    MEL(ndof_, W, 2, 1) = z_*(k_[2]*2/3*sqrt2*sy - kappa_[2]*2/3*sqrt2*sy)*cI;
#endif
  }
    
  /* This is a pair potential and so V0=V1=V  W0=W1=W */
  memcpy(V0, V, 9*sizeof(double_complex)); // V0

  if (nu_ > 1) {
    memcpy(V1, V, 9*sizeof(double_complex));
    memcpy(W0, W, 9*sizeof(double_complex));
    memcpy(W1, W, 9*sizeof(double_complex));
  }

  // Now need a U2, V2, W2, X0, X1, X2, X
  if (nu_ > 2) {
    /* nu_=1,2,4 supported, but interaction length only 1 or 2 supported.  
       Therefore U2=U3=U  V2=V3=V  W2=W3=W.  if nu_ > 2 copy bulk values */
    memcpy(U2, U, 9*sizeof(double_complex));
    memcpy(V2, V, 9*sizeof(double_complex));
    memcpy(W2, W, 9*sizeof(double_complex));
  }

  // Now need a U3, V3, W3, X3, Y0, Y1, Y2, Y3, Y
  if (nu_ > 3) {
    memcpy(U3, U, 9*sizeof(double_complex));
    memcpy(V3, V, 9*sizeof(double_complex));
    memcpy(W3, W, 9*sizeof(double_complex));
  }

  if (nu_ == 2 || nu_ == 4) {
    xcshift[1]=0.5; ycshift[1]=0.5; // fac[1] = cexp(COMPLEX_NUMBER(0.0, -mulphase_*(qx+qy)));
    xcshift[0]=0.0; ycshift[0]=0.0; // fac[0] = cexp(COMPLEX_NUMBER(0.0, mulphase_*(qx+qy)));
  }  

  if (nu_ == 4) {
    xcshift[3]=0.5; ycshift[3]=0.5; //  fac[3] = cexp(COMPLEX_NUMBER(0.0, -mulphase_*(qx+qy)));
    xcshift[2]=0.0; ycshift[2]=0.0; //  fac[2] = cexp(COMPLEX_NUMBER(0.0, mulphase_*(qx+qy)));
  }  

}

void FCC100StiffnessKernel::get_force_at_gamma_point(double *f)
{
  for (int i = 0; i < nu_; i++) {
    f[i] = linf_[i];
  }
}

/* ----------------------------------------------------------------------
 * Face centered cubic (100) surface, two atoms per unit cell
 * --------------------------------------------------------------------*/

FCC100x2StiffnessKernel::FCC100x2StiffnessKernel(int narg, int *carg,
						 char **arg,
						 Domain *domain,
						 Memory *memory,
						 Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  strcpy(name_, "fcc100x2");

  if (*carg >= narg) {
    error_->all(FLERR,"FCC100x2StiffnessKernel::FCC100x2StiffnessKernel: "
		"Expected number of interacting neighbors.");
  }

  nu_ = strtol(arg[*carg], &endptr, 10);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"FCC100x2StiffnessKernel::FCC100x2StiffnessKernel: Error "
		"converting dimension argument to integer.");
  }
  (*carg)++;

  if (nu_ != 1) {
    error_->all(FLERR,"FCC100x2StiffnessKernel::FCC100x2StiffnessKernel: "
		"Number of layers must be 1 for now.");
  }

  for (int i = 0; i < MIN(nu_, 2); i++) {
    if (*carg >= narg) {
      error_->all(FLERR,"FCC100x2StiffnessKernel::FCC100x2StiffnessKernel: "
		  "Spring constant expected.");
    }

    k_[i] = strtod(arg[*carg], &endptr);
    if (endptr == arg[*carg]) {
      error_->all(FLERR,"FCC100x2StiffnessKernel::FCC100x2StiffnessKernel: "
		  "Error converting spring constant to number.");
    }
    (*carg)++;
  }

  nu_ = 1;
  ndof_ = 6;
  dim_ = ndof_*nu_;
}

FCC100x2StiffnessKernel::~FCC100x2StiffnessKernel()
{
}

void
FCC100x2StiffnessKernel::get_per_layer_dynamical_matrices(double qx,
							  double qy,
							  double_complex **D,
							  double *ycshift, double *xcshift)
{
  double cx2, cy2, sx2, sy2;
  double_complex U0xx[3][3], Uxx[3][3], Uxy[3][3], Vxx[3][3], Vxy[3][3];
  double_complex *kU0, *kU, *kV;

  cx2 = cos(qx/2);
  cy2 = cos(qy/2);
  sx2 = sin(qx/2);
  sy2 = sin(qy/2);

  if (nu_ != 1) {
    error_->all(FLERR,"FCC100x2StiffnessKernel: nu_ != 1\n");
  }
  if (ndof_ != 6) {
    error_->all(FLERR,"FCC100x2StiffnessKernel: ndof_ != 6\n");
  }

  memset(U0xx, 0, 9*sizeof(double_complex));
  memset(Uxx,  0, 9*sizeof(double_complex));
  memset(Uxy,  0, 9*sizeof(double_complex));
  memset(Vxx,  0, 9*sizeof(double_complex));
  memset(Vxy,  0, 9*sizeof(double_complex));

  /* Uxx is diagonal */
  Uxx[0][0]  = k_[0]*4;
  Uxx[1][1]  = k_[0]*4;
  Uxx[2][2]  = k_[0]*4;

  /* U0xx is diagonal */
  U0xx[0][0] = Uxx[0][0] - k_[0]*1;
  U0xx[1][1] = Uxx[1][1] - k_[0]*1;
  U0xx[2][2] = Uxx[2][2] - k_[0]*2;

  /* diagonal of Uxy */
  Uxy[0][0]  = -k_[0]*2*cx2*cy2;
  Uxy[1][1]  = -k_[0]*2*cx2*cy2;

  /* off-diagonal of Uxy */
  Uxy[0][1]  =  k_[0]*2*sx2*sy2;
  Uxy[1][0]  =  k_[0]*2*sx2*sy2;

  /* diagonal of Vxx */
  Vxx[1][1]  = -k_[0]*cy2;
  Vxx[2][2]  = -k_[0]*cy2;

  /* off-diagonal of Vxx */
  Vxx[1][2]  =  k_[0]*sy2*cI;
  Vxx[2][1]  =  k_[0]*sy2*cI;

  /* diagonal of Vxy */
  Vxy[0][0]  = -k_[0]*cx2;
  Vxy[2][2]  = -k_[0]*cx2;

  /* off-diagonal of Vxy */
  Vxy[0][2]  =  k_[0]*sx2*cI;
  Vxy[2][0]  =  k_[0]*sx2*cI;

  /*
   * Assemble into kernel matrices
   */
  kU0 = D[0];
  kU  = D[1];
  kV  = D[2];

  memset(kU0, 0, dim_*dim_*sizeof(double_complex));
  memset(kU,  0, dim_*dim_*sizeof(double_complex));
  memset(kV,  0, dim_*dim_*sizeof(double_complex));

  put_matrix     (6, kU0, 3, &U0xx[0][0], 0, 0);
  put_matrix     (6, kU0, 3, &U0xx[0][0], 3, 3);
  put_matrix     (6, kU0, 3, &Uxy[0][0],  0, 3);
  put_matrix     (6, kU0, 3, &Uxy[0][0],  3, 0);

  put_matrix     (6, kU,  3, &Uxx[0][0],  0, 0);
  put_matrix     (6, kU,  3, &Uxx[0][0],  3, 3);
  put_matrix     (6, kU,  3, &Uxy[0][0],  0, 3);
  put_matrix     (6, kU,  3, &Uxy[0][0],  3, 0);

  put_matrix     (6, kV,  3, &Vxx[0][0],  0, 0);
  put_matrix     (6, kV,  3, &Vxx[0][0],  3, 3);
  put_matrix     (6, kV,  3, &Vxy[0][0],  0, 3);
  put_matrix     (6, kV,  3, &Vxy[0][0],  3, 0);
}

} /* namespace */
