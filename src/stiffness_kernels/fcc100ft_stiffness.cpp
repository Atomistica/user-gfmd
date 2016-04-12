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
// I see that fcc100ft is the same to 1% only as fcc100

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "fcc100ft_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Face centered cubic (100) surface
 * --------------------------------------------------------------------*/

FCC100FTStiffnessKernel::FCC100FTStiffnessKernel(int narg, int *carg, char **arg,
					     Domain *domain, Memory *memory,
					     Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  strcpy(name_, "fcc100ft");

  z_ = 1;
  mulphase_ = 0.5;

  k_[0] = 0.0;
  k_[1] = 0.0;
  k_[2] = 0.0;
  kappa_[0] = 0.0;
  kappa_[1] = 0.0;
  kappa_[2] = 0.0;

  if (*carg >= narg) {
    error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Expected "
		"number of interacting neighbors.");
  }

  // Hijacking the nu_ spot to find optional key word 'x' (ending in x)
  // Choosing to set kappa values
  int l = strlen(arg[*carg]);
  bool x = false;
  // x denotes kappa to be set 
  if (arg[*carg][l-1] == 'x') {
    x = true;
    printf("Setting the x option for kappa\n");
    (*carg)++;
    //arg[*carg][l] = '\0';
  }

  // Number of layers: nu_
  nu_ = strtol(arg[*carg], &endptr, 10); 
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Error "
		"converting number of layers to integer.");
  }
  (*carg)++;

  if (nu_ < 1)
    error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Number "
		"of bounday GF layers must be at least 1.");

  // Number of spring constants defined for the interacting neighbors: nk_
  switch( nu_ )
  {
    case 1: nk_ = 1; break; // 1 layer supports up to nn in the fcc lattice
    case 2: nk_ = 3; break; // 2 layers supports nn, 2n, and 3n spring constants
    default: nk_ = 3; 
  }

  for (int i = 0; i < nk_; i++) {
    if (*carg >= narg) {
      error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Spring "
		  "constant expected (1 for 1 layer, 3 for 2 layers.  Can be 0.");
    }

    k_[i] = strtod(arg[*carg], &endptr);
    if (endptr == arg[*carg]) {
      char errstr[1024];
      sprintf(errstr, "FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Expecting "
                      " %d ks.  Error converting spring constant to number.", nk_);
      error_->all(FLERR, errstr);
    }
    (*carg)++;

    if (x) {
      if (*carg >= narg) {
	error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: "
		    "Spring constant expected.");
      }

      kappa_[i] = strtod(arg[*carg], &endptr);
      if (endptr == arg[*carg]) {
        char errstr[1024];
        sprintf(errstr, "FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Expecting "
                      " %d ks.  Error converting spring constant to number.", nk_);
        error_->all(FLERR, errstr);
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
      height_ = strtol(arg[*carg], &endptr, 10);
      if (endptr == arg[*carg]) {
	error_->all(FLERR,"FCC100FTStiffnessKernel::FCC100FTStiffnessKernel: Error "
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
    fprintf(logfile, "FCC100FT boundary with %i layers and normal in direction "
	    "%i.", nu_, z_);
  }
  */

  // natomsunitcell_ = nu_; // For now assume one atom per layer in the unit cell
  dim_ = ndof_*nu_; // Tot dim=3*nc (in multiatom cell) or 3*num (1-atom cell)
  npars_ = 2*nk_; // Number of parameters, obsolete
}

FCC100FTStiffnessKernel::~FCC100FTStiffnessKernel()
{
  delete [] linf_;
}

void FCC100FTStiffnessKernel::set_parameters(double *pars)
{
  for (int i = 0; i < nk_; i++) {
    k_[i] = pars[i];
    kappa_[i] = pars[nk_+i];
  }
}

void FCC100FTStiffnessKernel::dump_info(FILE *f)
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


int FCC100FTStiffnessKernel::disttoconnecttype(double dist2) {
  // these dont need to be double but then need to be cast in later equations anyway
  double bordermax = 1.95/1.1; // This should come from somewhere that also gives is to fcc100fteam get_per...

  double border2to3 = 1.6/1.1;
  double border1to2 = 1.3/1.1;
  double bordermin = 1.0/1.1;
  if (dist2 > bordermax * bordermax) return -1;
  else if (dist2 > border2to3 * border2to3) return 2;
  else if (dist2 > border1to2 * border1to2) return 1;
  else if (dist2 > bordermin * bordermin ) return 0;
  return -1;
}


inline void setvec3(double *tar, double x, double y, double z)
{
  tar[0] = x;
  tar[1] = y;
  tar[2] = z;
}


// Multi-atom unit cells written as single atom cells, per Appendix A.3
// Will return a U V and W at this q vector
void FCC100FTStiffnessKernel::get_per_layer_dynamical_matrices(double qx,
							     double qy,
							     double_complex **D,
							       double *xcshift, double *ycshift)
{

#define NUMLNS 42 // Number lattice neighbors: 12 nn + 6 2nn + 24 3nn

  double sqrt2 = sqrt(2.0);
  double invsqrt2 = 1.0 / sqrt2;
  // wavevector (qx,qy) [2 pi / (lambda/spacing)]
  // grid spacing is also nearest neighbor distance = units of r below

  // D: list of pointers to the diagonal (ie non-zero) dynam matrix elements
  // If there are nu GF boundary layers, consider them as nu single-atom cells
  // D = { U0,(U1,...,Unu-1),U,V,(W,..) }
  // To match paper notation, nu-1 tildas sit above each variable name U, V...

  double_complex *U0, *U, *V, *U1, *W;
  // Handle up to nu = 4, but only interaction distance =2 so 
  // these are zero or bulk values:
  double_complex *U2, *U3, *V0, *V1, *V2, *V3, *W0, *W1, *W2, *W3;
  // D, to be general enough for multibody potentials will store all this:
  // U0 U1 U  (for nu_ == 2)   U0 U (for nu_ == 1)
  // V0 V1 V                   V0 V
  // W0 W1 W           
  //printf("nu_ %d", nu_);
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


  // Even if nu_ is four, only 2nd-nearest neighbors are included in the stiffness.
  // memset(D[0], 0, (2*nu_+1)*9*sizeof(double_complex));
  memset(D[0], 0, (nu_+1)*(nu_+1)*9*sizeof(double_complex));

  double para_k;
  double perp_k;
  double_complex phase = 0.0;
  SquareMatrix<double_complex> eij_eij(3);
  SquareMatrix<double_complex> perpeij_eij(3);
  SquareMatrix<double_complex> Umat(3);
  SquareMatrix<double_complex> Vmat(3);
  SquareMatrix<double_complex> Wmat(3);
  SquareMatrix<double_complex> contrib(3);
  SquareMatrix<double_complex> ucontribfromablow(3);
  SquareMatrix<double_complex> ucontribfromab1(3);
  SquareMatrix<double_complex> ucontribfromab2(3);

  // Since all data gets accumulated by complex matricies anyhow, use complex everywhere
  double_complex ident_data[9];
  for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
      ident_data[ii*3 + jj] = (ii == jj);
  SquareMatrix<double_complex> ident(3,ident_data);

  struct latticeneighbor {
    double r[3]; // Neighbors relative position [2pi/q] (since grid spacing = nn dist)
    double rnorm; // Precomputed length of r
    int kt; // Spring constant type. Convention: 0 for nearest neighbor type, 1 for 2nn...
    int ab; // alpha minus beta: relative z-position of this neighbor, units of atomic planes
  } ln [NUMLNS];
  double rij[3];

  // ab = 0
  // Easy-to-read "compound literal" syntax was allowed starting with C99 standard
  setvec3(ln[0 ].r,  1, 0,0);
  setvec3(ln[1 ].r,  0, 1,0);
  setvec3(ln[2 ].r, -1, 0,0);
  setvec3(ln[3 ].r,  0,-1,0);
  setvec3(ln[4 ].r,  1, 1,0);
  setvec3(ln[5 ].r,  1,-1,0);
  setvec3(ln[6 ].r, -1, 1,0);
  setvec3(ln[7 ].r, -1,-1,0);
  // ab = 1
  setvec3(ln[8 ].r,  0.5, 0.5, invsqrt2);
  setvec3(ln[9 ].r, -0.5, 0.5, invsqrt2);
  setvec3(ln[10].r,  0.5,-0.5, invsqrt2);
  setvec3(ln[11].r, -0.5,-0.5, invsqrt2);
  setvec3(ln[12].r,  1.5, 0.5, invsqrt2);
  setvec3(ln[13].r,  0.5, 1.5, invsqrt2);
  setvec3(ln[14].r, -1.5, 0.5, invsqrt2);
  setvec3(ln[15].r, -0.5, 1.5, invsqrt2);
  setvec3(ln[16].r,  1.5,-0.5, invsqrt2);
  setvec3(ln[17].r,  0.5,-1.5, invsqrt2);
  setvec3(ln[18].r, -1.5,-0.5, invsqrt2);
  setvec3(ln[19].r, -0.5,-1.5, invsqrt2);
  // ab = -1
  setvec3(ln[20].r,  0.5, 0.5,-invsqrt2);
  setvec3(ln[21].r, -0.5, 0.5,-invsqrt2);
  setvec3(ln[22].r,  0.5,-0.5,-invsqrt2);
  setvec3(ln[23].r, -0.5,-0.5,-invsqrt2);
  setvec3(ln[24].r,  1.5, 0.5,-invsqrt2);
  setvec3(ln[25].r,  0.5, 1.5,-invsqrt2);
  setvec3(ln[26].r, -1.5, 0.5,-invsqrt2);
  setvec3(ln[27].r, -0.5, 1.5,-invsqrt2);
  setvec3(ln[28].r,  1.5,-0.5,-invsqrt2);
  setvec3(ln[29].r,  0.5,-1.5,-invsqrt2);
  setvec3(ln[30].r, -1.5,-0.5,-invsqrt2);
  setvec3(ln[31].r, -0.5,-1.5,-invsqrt2);
  // ab = 2
  setvec3(ln[32].r,  0, 0, sqrt2);
  setvec3(ln[33].r,  1, 0, sqrt2);
  setvec3(ln[34].r,  0, 1, sqrt2);
  setvec3(ln[35].r, -1, 0, sqrt2);
  setvec3(ln[36].r,  0,-1, sqrt2);
  // ab = -2
  setvec3(ln[37].r,  0, 0,-sqrt2);
  setvec3(ln[38].r,  1, 0,-sqrt2);
  setvec3(ln[39].r,  0, 1,-sqrt2);
  setvec3(ln[40].r, -1, 0,-sqrt2);
  setvec3(ln[41].r,  0,-1,-sqrt2);

  // Distance between alpha=n and alpha=n+1 (single atom unit cells)
  double layerzseparation = invsqrt2; // should come from somewhere that also gives is to fcc100fteam get_per...

  // Useful to have rnorm precomputed.  
  // And need disttoconnecttype() for rjk anyway.
  for (int ii = 0; ii < NUMLNS; ii++) {
    ln[ii].rnorm = sqrt(vec_norm2(3,ln[ii].r));
    ln[ii].kt    = disttoconnecttype(ln[ii].rnorm * ln[ii].rnorm);
    ln[ii].ab    = round(ln[ii].r[2] / layerzseparation);
  }


  //  if ((nu_ > 1) && (fabs(k_[1]) > 1e-12) && (fabs(linf_[0]) < 1e-12)) {
  //    printf("Warning: linf_ not set despite greater than nearest neighbor"
  //           " interactions and so surface relaxation.\n");
  //    printf("Will set automatically, as can be done in fcc100ft if dphidr provided.\n");
  //    // Calculate force "matrix" which must be only q=0 and in the vertical direction only
  //    for (int kk = 0; kk < NUMLNS; kk++) {   
  //       // Pair potential force // Notes 4/28
  //      linf_[0] -=  dphidr_[ln[kk].kt] * (ln[kk].r)[2] / ln[kk].rnorm * (ln[kk].ab <=   0);
  //      if (nu_ > 1) linf_[1] -=  dphidr_[ln[kk].kt] * (ln[kk].r)[2] / ln[kk].rnorm * (ln[kk].ab <=   1);
  //    }  
  //  }

  // printmat(3,Umat.data()); // TAS3
  // printmat(3,Vmat.data()); // TAS3
  // printmat(3,Wmat.data()); // TAS3
  // Caution: kappa may need a minus sign depending on if it is PRB (A6) or -(A6)
  // Perform sums of Dij over atoms j (from the Fourier transform)
  // jj loops over all neighbors of i
  // Fix atom i since all Dij(r) are equivalent (except by distance from surface)
  for (int jj = 0; jj < NUMLNS; jj++) {

    setvec3(rij, -(ln[jj].r)[0],-(ln[jj].r)[1],-(ln[jj].r)[2]); // negative of rji 

    unit_outer_complex(3,rij,rij,eij_eij.data());
    perpeij_eij = ident - eij_eij;
    para_k = k_[ln[jj].kt];
    perp_k = kappa_[ln[jj].kt];

    // Assuming first two components of rel pos vector gives unit cell translation
    // Could extend operators available in linearalgebra.h or break down this expression:
    // contrib = (para_k * eij_eij + perp_k * perpeij_eij ) * exp(-I *( qx * rij[0] + qy * rij[1] ));
    phase = -1.0 * (qx * rij[0] + qy * rij[1]) * I;
    contrib = ( (eij_eij * (double_complex)(-para_k)) 
	+ (perpeij_eij * (double_complex)(-perp_k))) * exp(phase); // TAS3 something here is going wrong that affects q=0 V and W like r->-r that changes the sum magnitude 10% of x-x 

    switch (ln[jj].ab) {
      case  0: Umat += contrib; break; 
      case -1: Vmat += contrib; break; 
      case -2: Wmat += contrib; break; 
    }
  }

  // Final term in the sum over j: atom j = atom i: the delta function term
  // eik in the paper = eij in the code so can reuse matrix eij_eij
  for (int kk = 0; kk < NUMLNS; kk++) {
    para_k = k_[ln[kk].kt];
    perp_k = kappa_[ln[kk].kt];
    setvec3(rij, -(ln[kk].r)[0],-(ln[kk].r)[1],-(ln[kk].r)[2]); // negative of rki
    unit_outer_complex(3,rij,rij,eij_eij.data());
    perpeij_eij = ident - eij_eij;

    if (ln[kk].ab <= 0) {
      ucontribfromablow += eij_eij * ((double_complex)(para_k));
      ucontribfromablow += perpeij_eij * ((double_complex)(perp_k));
    } else if (ln[kk].ab == 1) {
      ucontribfromab1 += eij_eij * ((double_complex)(para_k));
      ucontribfromab1 += perpeij_eij * ((double_complex)(perp_k));
    } else if (ln[kk].ab == 2) {
      ucontribfromab2 += eij_eij * ((double_complex)(para_k));
      ucontribfromab2 += perpeij_eij * ((double_complex)(perp_k));
    }
  }

  // Copy the results into the D array pointers (U0, U, V...)
  memcpy(U0, (Umat + ucontribfromablow).data(), 9*sizeof(double_complex));
  memcpy(U, (Umat + ucontribfromablow + ucontribfromab1 + ucontribfromab2).data(), 9*sizeof(double_complex));
  memcpy(V, Vmat.data(), 9*sizeof(double_complex));
  memcpy(V0,Vmat.data(), 9*sizeof(double_complex));

  // For pair potentials V0 = V1 = V, and same with W, X, ...
  // Now need a U1, V1, W0, W1, W
  if (nu_ > 1) {
    memcpy(U1, (Umat + ucontribfromablow + ucontribfromab1).data(), 9*sizeof(double_complex));
    memcpy(W, Wmat.data(), 9*sizeof(double_complex));
    memcpy(V1, V, 9*sizeof(double_complex));
    memcpy(W0, W, 9*sizeof(double_complex));
    memcpy(W1, W, 9*sizeof(double_complex));
  }
  // Now need a U2, V2, W2, X0, X1, X2, X
  if (nu_ > 2) {
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

  // Set up for nu = 1, 2, 4 only (number of GF-tracked layers)
  if (nu_ == 2 || nu_ == 4) {
    xcshift[1]=0.5; ycshift[1]=0.5; // fac[1] = cexp(COMPLEX_NUMBER(0.0, -mulphase_*(qx+qy)));
    xcshift[0]=0.0; ycshift[0]=0.0; // fac[0] = cexp(COMPLEX_NUMBER(0.0, mulphase_*(qx+qy)));
  }  

  if (nu_ == 4) {
    xcshift[3]=0.5; ycshift[3]=0.5; //  fac[3] = cexp(COMPLEX_NUMBER(0.0, -mulphase_*(qx+qy)));
    xcshift[2]=0.0; ycshift[2]=0.0; //  fac[2] = cexp(COMPLEX_NUMBER(0.0, mulphase_*(qx+qy)));
  }  

}

void FCC100FTStiffnessKernel::get_force_at_gamma_point(double *f)
{
  for (int i = 0; i < nu_; i++) {
    f[i] = linf_[i];
  }
}


} /* namespace */


