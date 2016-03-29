#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "sc100_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
 * Simple cubic (100) surface, direct computation of the stiffness matrix
 * See: Y. Saito, J. Phys. Soc. Jpn. 73, 1816 (2004)
 * --------------------------------------------------------------------*/

SC100ExplicitStiffnessKernel::SC100ExplicitStiffnessKernel(int narg,
							   int *carg,
							   char **arg,
							   Domain *domain, 
							   Memory *memory,
							   Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  dim_ = 3;
  strcpy(name_, "sc100_explicit");
}

SC100ExplicitStiffnessKernel::~SC100ExplicitStiffnessKernel()
{
}

void SC100ExplicitStiffnessKernel::get_stiffness_matrix(double qx, double qy,
							double_complex *phi)
{
  double_complex U_q[3][3];
  double_complex F_q[3][3];

  if (qx == 0) qx = 1.e-6;
  if (qy == 0) qy = 1.e-6;
  if ( (fabs(qx) == M_PI) && (fabs(qy) == M_PI) )  qx=qy=M_PI-1.e-6;

  double c_x = cos(qx);
  double c_y = cos(qy);

  double s_x = sin(qx);
  double s_y = sin(qy);

  double c_plus = 1.0+ c_x + c_y - 2.*c_x*c_y;
  double c_minus = 4.0*sin(qx/2.0)*sin(qy/2.0)*sqrt(1.0-c_x*c_y);

  double s_plus = sqrt(  ( sqrt(pow ( (pow(c_plus,2.0)-pow(c_minus,2.0)-1.),2.0) + 4.*pow(c_plus,2.0)*pow(c_minus,2.0) )
            + ( pow(c_plus,2.0)-pow(c_minus,2.0)-1.) )  / 2.0 )*SIGN(c_plus);
  double s_minus = sqrt( ( sqrt( pow((pow(c_plus,2.0)-pow(c_minus,2.0)-1.),2.0) + 4.*pow(c_plus,2.0)*pow(c_minus,2.0) )
            - ( pow(c_plus,2.0)-pow(c_minus,2.0)-1.) ) /2.0 )*SIGN(c_minus);

  double e_plus = c_plus - s_plus;
  double e_minus = c_minus - s_minus;
  double e_one = (c_x + c_y +c_x*c_y)/(4.-c_x*c_y+sqrt((4.+c_x+c_y)*( 4.-c_x-c_y-2.*c_x*c_y)));

  U_q[0][0] = COMPLEX_NUMBER(sin(qx/2.0)*cos(qy/2.0)/(c_x+c_x*c_y-2.) ,0.);
  U_q[1][0] = COMPLEX_NUMBER(-1.*cos(qx/2.0)*sin(qy/2.0)/(c_y + c_x*c_y - 2.),0.);
  U_q[2][0] = COMPLEX_NUMBER(cos(qx/2.0)*cos(qy/2.0)*(c_x-c_y)/((c_y + c_x*c_y - 2.)*(c_x + c_x*c_y - 2.)),0.);

  U_q[2][0] = COMPLEX_NUMBER(creal(U_q[2][0])*(sqrt((4.-c_x-c_y-2.*c_x*c_y)/(4. + c_x + c_y))), cimag(U_q[2][0]));

  U_q[0][1] = COMPLEX_NUMBER(s_x/(2.*(1.-c_x+pow(c_x,2.0)-c_x*c_y)),0.);
  U_q[1][1] = COMPLEX_NUMBER(s_y/(2.*(1.-c_y+pow(c_y,2.0)-c_x*c_y)),0.);
  U_q[2][1] = COMPLEX_NUMBER(  ( (2.-2.*c_x*c_y*c_plus-pow(c_minus,2.0))*s_plus + (c_minus*(c_plus-2.*c_x*c_y))*s_minus)
            /( pow((2.-2.*c_x*c_y*c_plus-pow(c_minus,2.0)),2.0) + pow(c_minus*(c_plus-2.*c_x*c_y),2.0)),0.);

  U_q[0][2] = COMPLEX_NUMBER(c_x*cos(qx/2.)*sin(qy/2.)/(sqrt(1.-c_x*c_y)*(1.-c_x+pow(c_x,2.)-c_x*c_y)), 0.);
  U_q[1][2] = COMPLEX_NUMBER(c_y*sin(qx/2.)*cos(qy/2.)/(sqrt(1.-c_x*c_y)*(1.-c_y+pow(c_y,2.)-c_x*c_y)), 0.);

  U_q[2][2] = COMPLEX_NUMBER(((2.-2.*c_x*c_y*c_plus-pow(c_minus,2.0))*s_minus -(c_minus*(c_plus-2.*c_x*c_y))*s_plus)/
           (pow((2.-2.*c_x*c_y*c_plus-pow(c_minus,2.0)),2.0) +pow((c_minus*(c_plus-2.*c_x*c_y)),2.0) ),0.);

  F_q[0][0] = COMPLEX_NUMBER((5.-2.*c_x-2.*c_x*c_y-e_one*c_x)*creal(U_q[0][0]) + 2*s_x*s_y*creal(U_q[1][0]) + s_x*e_one*creal(U_q[2][0]),0.);
  F_q[0][1] = COMPLEX_NUMBER((5.-2.*c_x-2.*c_x*c_y-e_plus*c_x)*creal(U_q[0][1]) + 2*s_x*s_y*creal(U_q[1][1]) + c_x*e_minus*creal(U_q[0][2])
           + s_x*(e_plus*creal(U_q[2][1])-e_minus*creal(U_q[2][2])),0.);
  F_q[0][2] = COMPLEX_NUMBER((5.-2.*c_x-2.*c_x*c_y-e_plus*c_x)*creal(U_q[0][2]) + 2*s_x*s_y*creal(U_q[1][2]) - c_x*e_minus*creal(U_q[0][1])
           + s_x*(e_minus*creal(U_q[2][1])+e_plus*creal(U_q[2][2])),0.);
  F_q[1][0] = COMPLEX_NUMBER(2.*s_x*s_y*creal(U_q[0][0]) + (5.-2.*c_y-2.*c_x*c_y-e_one*c_y)*creal(U_q[1][0]) + s_y*e_one*creal(U_q[2][0]),0.);
  F_q[1][1] = COMPLEX_NUMBER(2.*s_x*s_y*creal(U_q[0][1]) + (5.-2.*c_y-2.*c_x*c_y-e_plus*c_y)*creal(U_q[1][1]) + e_minus*c_y*creal(U_q[1][2])
           +  s_y*(e_plus*creal(U_q[2][1])-e_minus*creal(U_q[2][2])),0.);
  F_q[1][2] = COMPLEX_NUMBER(2.*s_x*s_y*creal(U_q[0][2]) + (5.-2.*c_y-2.*c_x*c_y-e_plus*c_y)*creal(U_q[1][2]) -
           e_minus*c_y*creal(U_q[1][1]) + s_y*(e_minus*creal(U_q[2][1])+e_plus*creal(U_q[2][2])),0.);
  F_q[2][0] = COMPLEX_NUMBER(-s_x*e_one*creal(U_q[0][0])-s_y*e_one*creal(U_q[1][0]) + creal(U_q[2][0])*(3.-(1.+c_x+c_y)*e_one),0.);

  F_q[2][1] = COMPLEX_NUMBER(-s_x*(e_plus*creal(U_q[0][1])-e_minus*creal(U_q[0][2])) - s_y*(e_plus*creal(U_q[1][1])-e_minus*creal(U_q[1][2]))
           + (3.-(1.+c_x+c_y)*e_plus)*creal(U_q[2][1]) + (1.+c_x+c_y)*e_minus*creal(U_q[2][2]),0.);

  F_q[2][2] = COMPLEX_NUMBER(-s_x*(e_minus*creal(U_q[0][1])+e_plus*creal(U_q[0][2]))-s_y*(e_plus*creal(U_q[1][2])+e_minus*creal(U_q[1][1])) +
           creal(U_q[2][2])*(3.-(1.+c_x+c_y)*e_plus) - e_minus*(1.+c_x+c_y)*creal(U_q[2][1]),0.);

  GaussJordan(3,F_q[0],error_);
  MatMulMat(3,U_q[0],F_q[0],phi);
  
  MEL(dim_, phi, 0, 2) = COMPLEX_NUMBER(0,      creal(MEL(dim_, phi, 0, 2)));
  MEL(dim_, phi, 2, 0) = COMPLEX_NUMBER(0, -1.0*creal(MEL(dim_, phi, 2, 0)));
  MEL(dim_, phi, 1, 2) = COMPLEX_NUMBER(0,      creal(MEL(dim_, phi, 1, 2)));
  MEL(dim_, phi, 2, 1) = COMPLEX_NUMBER(0, -1.0*creal(MEL(dim_, phi, 2, 1)));

  GaussJordan(3,phi,error_);
}



/* ----------------------------------------------------------------------
 * Simple cubic (100) surface, transfer matrix implementation
 * --------------------------------------------------------------------*/

SC100StiffnessKernel::SC100StiffnessKernel(int narg, int *carg, char **arg,
					   Domain *domain, Memory *memory,
					   Error *error)
  : StiffnessKernel(narg, carg, arg, domain, memory, error)
{
  char *endptr;

  dim_ = 3;
  strcpy(name_, "sc100");

  if (*carg < narg-1) {
    if (!strcmp(arg[*carg], "height")) {
      (*carg)++;
      height_ = strtol(arg[*carg], &endptr, 10);
      if (endptr == arg[*carg]) {
	error_->all(FLERR,"SC100StiffnessKernel::SC100StiffnessKernel: Error "
		    "converting height_ argument to number.");
      }
      (*carg)++;
    }
  }
}

SC100StiffnessKernel::~SC100StiffnessKernel()
{
}

void SC100StiffnessKernel::get_per_layer_dynamical_matrices(double qx,
							    double qy,
							    double_complex **D,
							    double_complex *fac)
{
  double cx, cy, sx, sy;
  double_complex *U0;
  double_complex *U;
  double_complex *V;

  if (nu_ != 1) {
    error_->all(FLERR,"SC100StiffnessKernel::get_per_layer_dynamical_matrices:"
		"Only works for nu == 1.");
  }

  U0 = D[0];
  U  = D[1];
  V  = D[2];

  cx = cos(qx);
  cy = cos(qy);
  sx = sin(qx);
  sy = sin(qy);

  memset(U,  0, 9*sizeof(double_complex));
  memset(U0, 0, 9*sizeof(double_complex));
  memset(V,  0, 9*sizeof(double_complex));

  /* diagonal of U */
  MEL(dim_, U, 0, 0) = COMPLEX_NUMBER(6-2*cx*(1+cy), 0);
  MEL(dim_, U, 1, 1) = COMPLEX_NUMBER(6-2*cy*(1+cx), 0);
  MEL(dim_, U, 2, 2) = COMPLEX_NUMBER(6, 0);

  /* off-diagonal of U */
  MEL(dim_, U, 0, 1) = COMPLEX_NUMBER(2*sx*sy, 0);
  MEL(dim_, U, 1, 0) = COMPLEX_NUMBER(2*sx*sy, 0);

  /* diagonal of U0 */
  MEL(dim_, U0, 0, 0) = COMPLEX_NUMBER(5-2*cx*(1+cy), 0);
  MEL(dim_, U0, 1, 1) = COMPLEX_NUMBER(5-2*cy*(1+cx), 0);
  MEL(dim_, U0, 2, 2) = COMPLEX_NUMBER(3, 0);

  /* off-diagonal of U0 */
  MEL(dim_, U0, 0, 1) = COMPLEX_NUMBER(2*sx*sy, 0);
  MEL(dim_, U0, 1, 0) = COMPLEX_NUMBER(2*sx*sy, 0);

  /* diagonal of V */
  MEL(dim_, V, 0, 0) = COMPLEX_NUMBER(-cx, 0);
  MEL(dim_, V, 1, 1) = COMPLEX_NUMBER(-cy, 0);
  MEL(dim_, V, 2, 2) = COMPLEX_NUMBER(-1-cx-cy, 0);

  /* off-diagonal of V */
  MEL(dim_, V, 0, 2) = COMPLEX_NUMBER(0, sx);
  MEL(dim_, V, 1, 2) = COMPLEX_NUMBER(0, sy);
  
  MEL(dim_, V, 2, 0) = COMPLEX_NUMBER(0, sx);
  MEL(dim_, V, 2, 1) = COMPLEX_NUMBER(0, sy);
}

} /* namespace */
