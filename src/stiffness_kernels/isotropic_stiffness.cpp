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
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "isotropic_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Isotropic elasticity
 * --------------------------------------------------------------------*/

IsotropicStiffnessKernel::IsotropicStiffnessKernel(int narg, int *carg,
						   char **arg,
						   Domain *domain,
						   Memory *memory,
						   Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  dim_ = 3;
  strcpy(name_, "isotropic");

  // These were the default values
  nu_ = 0.0;
  myu_ = 1.0;
  gamma_ = 1.0;

  if (*carg >= narg+2) {
    error_->all(FLERR,"IsotropicStiffnessKernel::IsotropicStiffnessKernel: "
		"Expected Poisson number, shear modulus and gamma parameters.");
  }

  nu_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't convert nu to float.");
  }
  (*carg)++;

  myu_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't convert myu to float.");
  }
  (*carg)++;

  gamma_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't  convert gamma to float.");
  }
  (*carg)++;

  if (gamma_ > 0.0 && gamma_ < 1e6) {
    height_ = 1;
  }
}

IsotropicStiffnessKernel::~IsotropicStiffnessKernel()
{
}

void IsotropicStiffnessKernel::get_stiffness_matrix(double qx, double qy,
                                                    double_complex *phi,
                                                    double_complex dU)
{
  if (dU != 0.0) {
    error_->all(FLERR,"dU != 0 not supported.");
  }

  double xprd = domain_->xprd;
  double yprd = domain_->yprd;

  if (yprd != xprd)
    error_->all(FLERR, "isotropic kernel work with square surfaces only.");

  // parameters pertaining to analytic isotropic Green's function expressions

  double alpha = 1.0/(2.0*(1-nu_));
  double beta = sqrt(pow(qx,2.0)+pow(qy,2.0));

  if ( (fabs(qx) <= 1.e-06) && (fabs(qy) <= 1.e-06) ) {
    if (invariant_) {
      memset(phi, 0, dim_*dim_*sizeof(double_complex));
    }
    else {
      // Is it ok that this is cast to real? (TAS)
      MEL(3, phi, 0, 0) = COMPLEX_NUMBER((alpha+0.5)*xprd/(gamma_*M_PI),
					 0.0);
      MEL(3, phi, 1, 1) = MEL(3, phi, 0, 0);
      MEL(3, phi, 2, 2) = COMPLEX_NUMBER(xprd/(gamma_*M_PI),
					 0.0);
      MEL(3, phi, 0, 1) = MEL(3, phi, 0, 2) = MEL(3, phi, 1, 0) = 
	MEL(3, phi, 1, 2) = MEL(3, phi, 2, 0) = MEL(3, phi, 2, 1) = 
	COMPLEX_NUMBER(0.0, 0.0);

      // This inverts the Greens Function to create Phi
      //      GaussJordan(3, phi, error);
      invert3x3(phi);
    }
  }
  else {
    MEL(3, phi, 0, 0) = 
      COMPLEX_NUMBER(1.0/(myu_*beta) - 
		     (pow(qx,2.0)*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0))),
		     0.0);
    MEL(3, phi, 1, 1) = 
      COMPLEX_NUMBER(1.0/(myu_*beta) -
		     (pow(qy,2.0)*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0))),
		     0.0);
    MEL(3, phi, 2, 2) =
      COMPLEX_NUMBER((1.0/(2.0*myu_*alpha*beta)),
		     0.0);
    MEL(3, phi, 0, 1) =
      COMPLEX_NUMBER(-qx*qy*(2.0-1.0/alpha)/(2.0*myu_*pow(beta,3.0)),
		     0.0);
    MEL(3, phi, 1, 0) = MEL(3, phi, 0, 1);
    MEL(3, phi, 0, 2) = 
      COMPLEX_NUMBER(0.0,
		     -qx*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0)));
    MEL(3, phi, 1, 2) = 
      COMPLEX_NUMBER(0.0,
		     -qy*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0)));
    MEL(3, phi, 2, 0) =
      COMPLEX_NUMBER(0.0,
		     qx*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0)));
    MEL(3, phi, 2, 1) =
      COMPLEX_NUMBER(0.0,
		     qy*(1.0-1.0/alpha)/(2.0*myu_*pow(beta,2.0)));

    // This inverts the Greens Function to create Phi
    //    GaussJordan(3, phi, error);
    invert3x3(phi);
  }
}



/* ----------------------------------------------------------------------
 * Isotropic elasticity, z degree of freedom only
 * --------------------------------------------------------------------*/

IsotropicZStiffnessKernel::IsotropicZStiffnessKernel(int narg, int *carg,
						   char **arg,
						   Domain *domain,
						   Memory *memory,
						     Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  dim_ = 1;
  strcpy(name_, "isotropic/z");

  // These were the default values
  nu_ = 0.0;
  myu_ = 1.0;
  gamma_ = 1.0;

  if (*carg >= narg+2) {
    error_->all(FLERR,"IsotropicZStiffnessKernel::IsotropicStiffnessKernel: "
		"Expected Poisson number, shear modulus and gamma parameters.");
  }

  nu_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicZStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't convert nu to float.");
  }
  (*carg)++;

  myu_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicZStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't convert myu to float.");
  }
  (*carg)++;

  gamma_ = strtod(arg[*carg], &endptr);
  if (endptr == arg[*carg]) {
    error_->all(FLERR,"IsotropicZStiffnessKernel::IsotropicStiffnessKernel: "
		"Can't convert gamma to float.");
  }
  (*carg)++;

  if (gamma_ > 0.0 && gamma_ < 1e6) {
    height_ = 1;
  }
}

IsotropicZStiffnessKernel::~IsotropicZStiffnessKernel()
{
}

void IsotropicZStiffnessKernel::get_stiffness_matrix(double qx, double qy,
						    double_complex *phi)
{
  double xprd = domain_->xprd;

  // parameters pertaining to analytic isotropic Green's function expressions

  double alpha = 1.0/(2.0*(1-nu_));
  double beta = sqrt(pow(qx,2.0)+pow(qy,2.0));

  if ( (fabs(qx) <= 1.e-06) && (fabs(qy) <= 1.e-06) ) {
    if (invariant_) {
      phi[0] = 0.0;
    }
    else {
      phi[0] = 1.0/COMPLEX_NUMBER(xprd/(gamma_*M_PI), 0.0);
    }
  }
  else {
    phi[0] = 1.0/COMPLEX_NUMBER((1.0/(2.0*myu_*alpha*beta)), 0.0);
  }
}

} /* namespace */
