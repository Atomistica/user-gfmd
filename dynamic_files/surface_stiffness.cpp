#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "surface_stiffness.h"

#include "debug_stiffness.h"
#include "isotropic_stiffness.h"

#ifdef GFMD_FFTW3
#include "nonperiodic_stiffness.h" 
#endif

#include "sc100_stiffness.h"
#include "fcc100_stiffness.h"
#include "fcc100ft_stiffness.h"

#ifdef GFMD_MANYBODY
#include "fcc100fteam_stiffness.h"
#endif

#include "ft_stiffness.h"

using namespace std;
using namespace LAMMPS_NS;

#define SIGN(a) ((a)>0.?1.:((a)<0.?-1.:0.))

namespace LAMMPS_NS {

double StiffnessKernel::unit_cell[9] = { 1,0,0, 0,1,0, 0,0,1 };

/* ----------------------------------------------------------------------
 * StiffnessKernel class
 * Computes dynamical matrices for each q-vector, and the associated
 * stiffness matrices
 * --------------------------------------------------------------------*/

StiffnessKernel::StiffnessKernel(int narg, int *carg, char **arg,
                                 Domain *domain, Memory *memory, 
                                 Error *error)
{
  domain_ = domain;
  memory_ = memory;
  error_ = error;

  ndof_ = 3;
  nu_ = 1;
  dim_ = ndof_*nu_;
  height_ = 0;
  npars_ = 0;
  invariant_ = false;
  strcpy(name_, "N/A");

  D_ = NULL;
}


StiffnessKernel::~StiffnessKernel()
{
  if (D_) {
    memory_->destroy(D_);
  }
}


/*!
 * Dump generic text info to file
 */
void StiffnessKernel::dump_info(FILE *f)
{
    fprintf(f, "Elasticity kernel = '%s', total substrate height = %i "
            "layers, degrees of freedom per layer = %i.\n", name_, height_,
            dim_);
}


void StiffnessKernel::set_parameters(double *)
{
  char errstr[80];
  sprintf(errstr, "Stiffness kernel '%s' does not provide a mechanism to "
	  "change parameters after instantiation.", name_);
  error_->all(FLERR,errstr); 
}


void StiffnessKernel::get_per_layer_dynamical_matrices(double qx, double qy,
						       double_complex **D,
						       double *xcshift, double *ycshift)
{
  //printf("\n Entering get per layer dynamical matrices \n"); // dbg_JMM
  char errstr[80];
  sprintf(errstr, "Stiffness kernel '%s' does not provide per-layer dynamical "
          " matrices.", name_);
  error_->all(FLERR,errstr); 
}


/*
 * Automatically assemble dynamical matrices for larger unit cells based on
 * the single unit cell ones obtained by get_per_layer_...
 */
void StiffnessKernel::get_dynamical_matrices(double qx, double qy,
                                             double_complex *U0,
                                             double_complex *U,
                                             double_complex *V,
                                             double_complex dU)
{
  //printf("\n Entering get dynamical matrices\n"); // dbg_JMM
  //printf("\n qx %f | qy %f \n",qx,qy); // dbg_JMM
  int dim_sq;
  int i_D, j_D;
  double_complex *m1, *m2, *m3;
  double xcshift[nu_];
  double ycshift[nu_];
  double_complex fac;
  //printf("\n dim %i ndof_ %i nu_ %i\n",dim_, ndof_,nu_); // dbg_JMM
  
  if (dim_ != ndof_*nu_) {
    error_->all(FLERR,"Internal: dim != ndof_*nu_ in "
		"StiffnessKernel::get_dynamical_matrices");
  }

  for (int i = 0; i < nu_; i++) {
    xcshift[i] = 0.0;
    ycshift[i] = 0.0;
  }

  if (!D_) {
    /* D is a list of matrices of dimension ndof_ x ndof_ */
    int nummatricies =  (nu_ + 1)*(nu_ + 1); // ((nu_*nu_) +(3*nu_) +2)/2;
    memory_->create(D_, nummatricies, ndof_*ndof_, "StiffnessKernel::D");
  }

  /*
   * D_ array holds...
   * ...for nu_ = 1 -> (2x2)     U0, U,
   *                             V0, V        with tildas implied
   * ...for nu_ = 2 -> (3x3)     U0, U1, U,   
   *                             V0, V1, V,  
   *                             W0, W1, W    with 2xtildas implied
   * ...for nu_     ->(nu_+1)^2 
   */
  
  get_per_layer_dynamical_matrices(qx, qy, D_, xcshift, ycshift);

  dim_sq = dim_*dim_;
  memset(U0, 0, dim_sq*sizeof(double_complex));
  memset(U,  0, dim_sq*sizeof(double_complex));
  memset(V,  0, dim_sq*sizeof(double_complex));


  /*
   * Cell-Level U0, U, V0, and V matricies are composed of tilda'd matricies
   * For nu_ = 3 and      U= U' V' W'   V= X' 0' 0'   U0=U0' V0' W0'
   * interaction up to X     V" U' V'      W' X' 0'      V0" U1' V1'
   * ' = tildax2             W" V" U'      V' W' X'      W0" V1" U2'
   * " = tildax2 transpose             
   */

  // Traverse cell-level matrix U (and U0) going atom-to-atom
  // Example: fcc100 (1nn,2nn,3nn) lattice has one atom per nu_, and ndof=3
  // Following layout of PRB Eq 10, i counts like alpha, j from left to right
  // j-i is distance above diagonal in U
  // Diagrams in notes at bottom
  for (int i = 0; i < nu_; i++) {
    for (int j  = i; j < nu_; j++) {
      // Map from location in cell-level matrix U to D_ array
      i_D = j - i;
      j_D = i;
      m1 = D_[i_D * (nu_+1) + j_D];
      m2 = D_[i_D * (nu_+1) + nu_];
      // Copy that tilda'd matrix (e.g. U0' or U' or W3' or X3') into U or U0
      if (i == j) {
        for (int k = 0; k < ndof_; k++) {
          for (int l = 0; l < ndof_; l++) {
    	    MEL(dim_, U0, k+i*ndof_, l+j*ndof_) = MEL(ndof_, m1, k, l); 
    	    MEL(dim_, U,  k+i*ndof_, l+j*ndof_) = MEL(ndof_, m2, k, l);
          }
        }      
      } else { // i != j
	// Add in also the conjugate transpose element, opposite the diagonal
        // Fac is the phase modification from unit cell translations from layer i to j
        fac=cexp(COMPLEX_NUMBER(0.0,(xcshift[j]-xcshift[i])*qx +(ycshift[j]-ycshift[i])*qy));
	for (int k = 0; k < ndof_; k++) {
          for (int l = 0; l < ndof_; l++) { 
    	    MEL(dim_, U0, k+i*ndof_, l+j*ndof_) = fac *MEL(ndof_, m1, k, l);
            MEL(dim_, U,  k+i*ndof_, l+j*ndof_) = fac *MEL(ndof_, m2, k, l);
      	    MEL(dim_, U0, k+j*ndof_, l+i*ndof_) = conj(fac *MEL(ndof_, m1, l, k));
	    MEL(dim_, U,  k+j*ndof_, l+i*ndof_) = conj(fac *MEL(ndof_, m2, l, k));
          }
        }
      }
    }
  }

  // Traverse cell-level matrix V
  // Following layout of PRB Eq 10, i counts like alpha, j from left to right
  // i-j is dist below diagonal in V matrix
  for (int i = 0; i < nu_; i++) {
    for (int j  = 0; j < i+1; j++) {
      // Map from location in cell-level matrix V to D_ array
      i_D = nu_ + j - i;
      j_D = nu_;
      m1 = D_[i_D * (nu_+1) + j_D];
      // dbg_JMM
      //printf("\n i_D %d | j_d %d | index %d \n ",i_D,j_D, i_D * (nu_+1) + j_D); 
      for (int k = 0; k < ndof_; k++) { 
        for (int l = 0; l < ndof_; l++) {
          fac = cexp(COMPLEX_NUMBER(0.0,(xcshift[j]-xcshift[i])*qx +(ycshift[j]-ycshift[i])*qy));
	  MEL(dim_, V, k+i*ndof_, l+j*ndof_) =fac  *(MEL(ndof_, m1, k, l));
          // printf("wrote to %d \n", &(MEL(dim_, V,  k+i*nu_, l+j))  );
        }
      }      
    }
  }

  /*
   * Add mass and damping terms to U0 and U
   */

  for (int k = 0; k < dim_; k++) {
    U0[k*(dim_+1)] += dU;
    U[k*(dim_+1)] += dU;
  }

  
#if 0
  //if ((qx > .2) && (qx < 2.2) && (qy > .2) && (qy < 2.2)) {
  //printf("qx = %f, qy = %f\n", qx, qy);
  //printf("Original matrices...\n");
  for (int i = 0; i < nu_; i++) {
    printf("U%i:\n", i);
    printmat(ndof_, D_[i]);
  }
  printf("\nU:\n");
  printmat(ndof_, D_[nu_]);

  for (int i = 0; i < nu_; i++) {
    printf("V%i:\n", i);
    printmat(ndof_, D_[nu_+1 + i]);
  }
  printf("\nV:\n");
  printmat(ndof_, D_[2*nu_+1]);

  if (nu_ > 1) {
    for (int i = 0; i < nu_; i++) {
      printf("W%i:\n", i);
      printmat(ndof_, D_[2*nu_+2 + i]);
    }
    printf("\nW:\n");
    printmat(ndof_, D_[3*nu_+2]);
  }
  if (nu_ > 2) {
    for (int i = 0; i < nu_; i++) {
      printf("X%i:  is el D_[%d] \n", i, 3*nu_+3 + i);
      printmat(ndof_, D_[3*nu_+3 + i]);
    }
    printf("\nX:\n");
    printmat(ndof_, D_[4*nu_+3]);
  } 
#endif
#if 0
  printf("\nRenormalized matrices...\n");
  printf("U0:\n");
  printmat(dim_, U0);
  printf("U:\n");
  printmat(dim_, U);
  printf("V\n");
  printmat(dim_, V);
  printf("\n");
#endif
}


void StiffnessKernel::get_stiffness_matrix(double qx, double qy,
                                           double_complex *phi,
                                           double_complex dU)
{
  // printf("\n Entering get stiffness matrix 1\n"); // dbg_JMM
  SquareMatrix<double_complex> U0(dim_), U(dim_), V(dim_);

  get_dynamical_matrices(qx, qy, U0.data(), U.data(), V.data(), dU);

  // dbg_JMM
  
#if 0
    printf("qx = %f, qy = %f\n", qx, qy);
    printf("U0 =\n");
    U0.print();
    printf("U =\n");
    U.print();
    printf("V =\n");
    V.print();
#endif

  // use faster renormalization group calculcation if height is power of 2
  // Note: This breaks if height_ = 0. Fix when re-enabling.
  // int M = lround(log(height_)/log(2.0));
  // int N = lround(pow(2.0, M));
  // TAS comment out this for initial testing; force use of transfer_matrix_stiffness
  //if (height_ == N) {
  //  renormalization_group_stiffness(height_, ndof_, dim_, U0.data(), U.data(),
  //				      V.data(), phi, error_);
  //}
  //else {
    //        printf(" \n phi \n");
    //   printmat(ndof_,phi); // dbg_JMM
  greens_function_transfer_matrix_stiffness(height_, ndof_, dim_, U0.data(),
                                            U.data(), V.data(), phi, error_);

  //}
 if (invariant_ && fabs(qx) < 1e-12 && fabs(qy) < 1e-12) {
    printf("TAS2 we're enforcing sum rule\n");
    enforce_phi0_sum_rule(ndof_, phi, 1e-15);
  }

  // TAS2 Does phi match what you expect? For this height parameter?
  //if((qx == 0.0) && (qy == 0.0))  {printf("TAS2 here is phi: \n"); printmat(dim_,phi);}// TAS2
  // --- DEBUG ---
#if 0
  printf("Transfer matrix:\n");
  printf("qx = %f, qy = %f\n", qx/(2*M_PI), qy/(2*M_PI));
  printf("height = %i\n", height_);
  printf("GF:\n");
  greens_function_transfer_matrix_stiffness(height_, ndof_, dim_, U0.data(),
					    U.data(), V.data(), phi, error_);
  printmat(dim_, phi);
  if (is_Hermitian(dim_, phi))
    printf("is Hermitian\n");
  else
    printf("is not Hermitian\n");
#endif

#if 0
  printf("qx = %f, qy = %f\n", qx/(2*M_PI), qy/(2*M_PI));
  printf("height = %i\n", height_);
  SquareMatrix<double_complex> phi2(dim_);
  renormalization_group_stiffness(height_, ndof_, dim_, U0.data(), U.data(),
				  V.data(), phi2.data(), error_);
  printf("Transfer matrix:\n");
  printmat(dim_, phi);
  printf("Renormalization group:\n");
  phi2.print();
#if 0
  if (is_Hermitian(dim_, phi))
    printf("is Hermitian\n");
  else
    printf("is not Hermitian\n");
#endif
  phi2 -= phi;
  if (phi2.norm2() > 1e-6) {
    error_->all(FLERR, "phis not identical.");
  }
#endif

#if 0
  printf("direct:\n");
  direct_inversion_stiffness(height_, ndof_, dim_, U0, U, V, phi, error_);
  printmat(dim_, phi);
  if (is_Hermitian(dim_, phi))
    printf("is Hermitian\n");
  else
    printf("is not Hermitian\n");
#endif

#if 0
  printf("displacements:\n");
  displacement_transfer_matrix_stiffness(height_, ndof_, dim_, U0, U, V, phi,
					 error_);
  printmat(dim_, phi);
  if (is_Hermitian(dim_, phi))
    printf("is Hermitian\n");
  else
    printf("is not Hermitian\n");

  exit(999);
#endif
  // --- END DEBUG ---

  
}

void StiffnessKernel::get_stiffness_matrix(int height, double qx, double qy,
					   double_complex *phi,
					   double_complex dU)
{
  // We want to truncate dynamic q-chains with q-dependent depth. Use height,
  // opposed to height_ when calculating phi
  SquareMatrix<double_complex> U0(dim_), U(dim_), V(dim_);
  
  get_dynamical_matrices(qx, qy, U0.data(), U.data(), V.data(), dU);
  
  greens_function_transfer_matrix_stiffness(height, ndof_, dim_, U0.data(),
					    U.data(), V.data(), phi, error_);
  
  if (invariant_ && fabs(qx) < 1e-12 && fabs(qy) < 1e-12) {
    printf("TAS2 we're enforcing sum rule\n");
    enforce_phi0_sum_rule(ndof_, phi, 1e-15);
  }
}

void StiffnessKernel::get_stiffness_matrix(int nx, int ny, double qx, double qy,
                                           double_complex *phi,
                                           double_complex dU)
{
  //  printf("\n Entering get stiffness matrix 2\n"); // dbg_JMM
  get_stiffness_matrix(qx, qy, phi, dU);
  // subclasses may implement get_stiffness_matrix(qx,...)
  //   (eg isotropic_stiffness does this).
  // Otherwise surface_stiffness provides default implementation 
  //   which calls get_dynamical w.c. get_per_layer_dynamical

}


void StiffnessKernel::get_Gn0(int n, double qx, double qy, double_complex *Gn0,
			      double_complex dU)
{
  SquareMatrix<double_complex> U0(dim_), U(dim_), V(dim_), phi(dim_);
  //if ((qx > 1 ) && (qx < 2) && (qy > 1) && (qy < 2)) printf("calling from get_Gn0\n"); // TAS
  get_dynamical_matrices(qx, qy, U0.data(), U.data(), V.data(), dU);
  greens_function_transfer_matrix_stiffness(height_, ndof_, dim_, U0.data(),
                                            U.data(), V.data(), phi.data(), error_);
  transfer_matrix_Gn0(n, height_, ndof_, dim_, U0.data(), U.data(), V.data(),
                      phi.data(), Gn0, error_);
}


void StiffnessKernel::get_force_at_gamma_point(double *f)
{
  for (int i = 0; i < nu_; i++) {
    f[i] = 0.0;
  }
}



/* ----------------------------------------------------------------------
 * Instantiate a stiffness kernel according to keyword arguments
 * --------------------------------------------------------------------*/

StiffnessKernel *stiffness_kernel_factory(char *keyword, int narg, int *carg,
                                          char **arg, Domain *domain,
                                          Force *force, Memory *memory,
                                          Error *error)
{
  //  printf("\n Entering stiffness kernel factory \n"); // dbg_JMM
  char errstr[120];
  StiffnessKernel *kernel = NULL;

  
  if (!strcmp(keyword, "debug")) {
    kernel = new DebugStiffnessKernel(narg, carg, arg, domain, memory, error);
  }
  else if (!strcmp(keyword, "chain")) {
    kernel = new ChainStiffnessKernel(narg, carg, arg, domain, memory, error);
  }
  else if (!strcmp(keyword, "isotropic")) {
    kernel = new IsotropicStiffnessKernel(narg, carg, arg, domain, memory,
					  error);
  }
  else if (!strcmp(keyword, "isotropic/z")) {
    kernel = new IsotropicZStiffnessKernel(narg, carg, arg, domain, memory,
					   error);
  }
#ifdef GFMD_FFTW3
  else if (!strcmp(keyword, "nonperiodic")) {
    kernel = new NonperiodicStiffnessKernel(narg, carg, arg, domain, memory,
                                            error);
  }
#endif
  else if (!strcmp(keyword, "sc100_explicit")) {

    printf("\n Now entering SC100 exp kernel"); // dbg_JMM
    kernel = new SC100ExplicitStiffnessKernel(narg, carg, arg, domain, memory,
					      error);
    printf("\n Now leaving surface_stiffness with kernel in tow \n"); // dbg_JMM
  }
  else if (!strcmp(keyword, "sc100")) {
    kernel = new SC100StiffnessKernel(narg, carg, arg, domain, memory, error);
  }
  else if (!strcmp(keyword, "fcc100")) {
    kernel = new FCC100StiffnessKernel(narg, carg, arg, domain, memory, error);
  }
  else if (!strcmp(keyword, "fcc100ft")) {
    kernel = new FCC100FTStiffnessKernel(narg, carg, arg, domain, memory, 
                                          error);
  }
#ifdef GFMD_MANYBODY
  else if (!strcmp(keyword, "fcc100fteam")) {
    kernel = new FCC100FTEAMStiffnessKernel(narg, carg, arg, domain, memory, 
                                          error);
  }
#endif
  else if (!strcmp(keyword, "ft")) {
      kernel = new FTStiffnessKernel(narg, carg, arg, domain, force, memory, error);
  }

  if (kernel) {
    if (strcmp(kernel->get_name(), keyword)) {
      sprintf(errstr, "stiffness_kernel_factory: Internal error: keyword '%s' "
	      "and kernel '%s' name mismatch.", keyword, kernel->get_name());
      error->all(FLERR,errstr);
    }
  }

  return kernel;
}



/* ----------------------------------------------------------------------
 * get the solution of the (quadratic) equation X = A1*(B + X)^-1*A2
 * --------------------------------------------------------------------*/
// GFMD calls this with:
//  X = VT_{N-1} where N is the specified finite height
//  B = U
//  A1= V
//  A2=-Vdagger 
void iterate_Gnn(int dim, double_complex *A1, double_complex *A2,
		 double_complex *B, double_complex *X, double eps,
		 int maxit, Error *error)
{
  //  printf(" \n Entering iterate Gnn \n"); // dbg_JMM
  int dim_sq = dim*dim;
#ifdef HAVE_C99
  double_complex last_X[dim_sq];
  double_complex tmp[dim_sq];
#else
  double_complex *last_X, *tmp;
  last_X = new double_complex[dim_sq];
  tmp = new double_complex[dim_sq];
#endif
  int it;
  double dnorm = 1.0;
  double eps2 = eps*eps;

  memcpy(last_X, X, dim_sq*sizeof(double_complex));

  it = 0;
  while (dnorm > eps2 && it < maxit) {

    /* X = X + B */
    for (int i = 0; i < dim_sq; i++) {
      X[i] += B[i];
    }

    /* X = X^-1 */
    GaussJordan(dim, X, error);

    /* tmp = X*A2 */
    MatMulMat(dim, X, A2, tmp);

    /* X = A1*tmp = A1*X*A2 */
    MatMulMat(dim, A1, tmp, X);

    /* check convergence */
    dnorm = 0.0;
    for (int i = 0; i < dim_sq; i++) {
      dnorm += cnorm(X[i] - last_X[i]);
    }
    memcpy(last_X, X, dim_sq*sizeof(double_complex));

    it++;

  }

  //  printf("\n Here's how many iterations that took: %d",it); // dbg_JMM
  
  if (it >= maxit && eps2 > 0.0) {
    error->all(FLERR,"Out of iterations while evaluating the continued fraction.");
  }

#ifndef HAVE_C99
  delete [] last_X;
  delete [] tmp;
#endif
}


#if 0
/* ----------------------------------------------------------------------
 * get the solution of the equation X = A1*(B + X)^-1*A2
 * --------------------------------------------------------------------*/
void iterate_Gn0(int dim, double_complex *V,
		 double_complex *Vdagger, double_complex *U,
		 double_complex *Gn1, double_complex *Gn2,
		 int maxit, Error *error)
{
  int dim_sq = dim*dim;
#ifdef HAVE_C99
  double_complex inv_V[dim_sq], tmp[dim_sq];
#else
  double_complex *inv_V, *tmp;
  inv_V = new double_complex[dim_sq];
  tmp = new double_complex[dim_sq];
#endif

  memcpy(inv_V, V, dim_sq*sizeof(double_complex));
  GaussJordan(dim, inv_V, error);

  for (int i = 0; i < dim_sq; i++) {
    inv_V[i] = -inv_V[i];
  }

  for (int i = 0; i < maxit; i++) {
    /* tmp = U*Gn1 */
    mat_mul_mat(dim, U, Gn1, tmp);
    /* tmp = tmp + Vdagger*Gn2 */
    mat_mul_mat(dim, 1.0, Vdagger, Gn2, 1.0, tmp);

    /* Gn2 = Gn1 */
    memcpy(Gn2, Gn1, dim_sq*sizeof(double_complex));
    /* Gn1 = -V^-1*tmp */
    mat_mul_mat(dim, inv_V, tmp, Gn1);
  }

#ifndef HAVE_C99
  delete [] inv_V;
  delete [] tmp;
#endif
}
#endif


// ndof is degrees of freedom per atom; 3 in 3D
void reverse_atom_order(int ndof, int dim, double_complex *out,
			double_complex *in)
{
  int nlayers = dim/ndof;

  for (int i = 0; i < nlayers; i++) {
    for (int j = 0; j < nlayers; j++) {
      for (int ii = 0; ii < ndof; ii++) {
	for (int jj = 0; jj < ndof; jj++) {
	  MEL(dim, out, (nlayers-i-1)*ndof+ii, (nlayers-j-1)*ndof+jj) = 
	    MEL(dim, in, i*ndof+ii, j*ndof+jj);
	}
      }
    }
  }
}



/* ----------------------------------------------------------------------
 * private method, to get the analytic stiffness coefficients for
 * using full matrix inversion
 * --------------------------------------------------------------------*/
void direct_inversion_stiffness(int height,
				int ndof,
				int dim,
				double_complex *U0,
				double_complex *U,
				double_complex *V,
				double_complex *phi,
				Error *error)
{
  int dim_sq = dim*dim;
  int full_dim = dim*(height+1);
  int full_dim_sq = full_dim*full_dim;
#ifdef HAVE_C99
  double_complex D[full_dim_sq], UN[dim_sq], G[4*dim_sq];
#else
  double_complex *D, *UN, *G;
  D = new double_complex[full_dim_sq];
  UN = new double_complex[dim_sq];
  G = new double_complex[4*dim_sq];
#endif

  if (height < 1)
    error->all(FLERR,"direct_inversion_stiffness: height < 1\n");

  //  reverse_atom_order(ndof, dim, UN, U0);
  memcpy(UN, U, dim_sq*sizeof(double_complex));

  memset(D, 0, full_dim_sq*sizeof(double_complex));

#if 0
  printf("U0\n");
  printmat(dim, U0);
  printf("UN\n");
  printmat(dim, UN);
  printf("U\n");
  printmat(dim, U);
  printf("V\n");
  printmat(dim, V);
#endif

  put_matrix(full_dim, D, dim, U0, 0, 0);
  put_matrix(full_dim, D, dim, UN, dim*height, dim*height);
  for (int i = 1; i < height; i++) {
    put_matrix     (full_dim, D, dim, U, dim*i,     dim*i    );
  }
  for (int i = 0; i < height; i++) { 
    put_matrix     (full_dim, D, dim, V, dim*i,     dim*(i+1));
    put_conj_matrix(full_dim, D, dim, V, dim*(i+1), dim*i    );
  }

  //  printf("D:\n");
  //  printmat(full_dim, D);
  GaussJordan(full_dim, D, error);
  //  printf("B\n");

#if 1
  take_matrix(full_dim, D, dim, phi, 0, 0);
  put_matrix(2*dim, G, dim, phi, 0, 0);
  take_matrix(full_dim, D, dim, phi, dim*height, dim*height);
  put_matrix(2*dim, G, dim, phi, dim, dim);
  take_matrix(full_dim, D, dim, phi, 0, dim*height);
  put_matrix(2*dim, G, dim, phi, 0, dim);
  take_matrix(full_dim, D, dim, phi, dim*height, 0);
  put_matrix(2*dim, G, dim, phi, dim, 0);

  printf("G (small):\n");
  printmat(2*dim, G);

  GaussJordan(2*dim, G, error);
  take_matrix(2*dim, G, dim, phi, 0, 0);
#endif

#if 0
  take_matrix(full_dim, D, dim, phi, 0, 0);
  GaussJordan(dim, phi, error);
#endif

#ifndef HAVE_C99
  delete [] D;
  delete [] UN;
  delete [] G;
#endif

  //  printf("C\n");
}


/* ----------------------------------------------------------------------
 * private method, to get the analytic stiffness coefficients for
 * using full matrix inversion
 * --------------------------------------------------------------------*/
void displacement_transfer_matrix_stiffness(int height,
					    int ndof,
					    int dim,
					    double_complex *U0,
					    double_complex *U,
					    double_complex *V,
					    double_complex *phi,
					    Error *error)
{
  int dim_sq = dim*dim;
  int dim2 = 2*dim;
  int dim2_sq = dim2*dim2;
#ifdef HAVE_C99
  double_complex inv_V[dim_sq], tmp1[dim_sq], tmp2[dim_sq];
  double_complex T[dim2_sq], TN[dim2_sq], tmp[dim2_sq];
#else
  double_complex *inv_V, *tmp1, *tmp2, *T, *TN, *tmp;
  inv_V = new double_complex[dim_sq];
  tmp1 = new double_complex[dim_sq];
  tmp2 = new double_complex[dim_sq];
  T = new double_complex[dim2_sq];
  TN = new double_complex[dim2_sq];
  tmp = new double_complex[dim2_sq];
#endif

  if (height < 1)
    error->all(FLERR,"displacement_transfer_matrix_stiffness: height < 1\n");

  //  printf("V:\n");
  //  printmat(dim, V);

  conj_transpose(dim, inv_V, V);
  //  printf("First GaussJordan\n"); 
  GaussJordan(dim, inv_V, error);
  //  printf("...okay\n"); 

  //  printf("inv_V:\n");
  //  printmat(dim, inv_V);

  mat_mul_mat(dim, inv_V, V, tmp1);
  //  printmat(dim, tmp1);

  memset(T, 0, dim2_sq*sizeof(double_complex));
  for (int i = 0; i < dim; i++)
    MEL(dim2, T, i, dim+i) = 1.0;

  mat_mul_mat(dim, -1.0, inv_V, V, 0.0, tmp1);
  put_matrix(dim2, T, dim, tmp1, dim, 0);
  mat_mul_mat(dim, -1.0, inv_V, U, 0.0, tmp1);
  put_matrix(dim2, T, dim, tmp1, dim, dim);

  //  printf("T:\n");
  //  printmat(dim2, T);

  /* TN = T^(N-1) */
  memcpy(TN, T, dim2_sq*sizeof(double_complex));
  for (int i = 0; i < height-1; i++) {
    mat_mul_mat(dim2, T, TN, tmp);
    memcpy(TN, tmp, dim2_sq*sizeof(double_complex));
  }

  //  printf("TN:\n");
  //  printmat(dim2, TN);

  take_matrix(dim2, TN, dim, tmp1, 0, dim);
  take_matrix(dim2, TN, dim, tmp2, dim, dim);

  //  printmat(dim, tmp1);
  //  printmat(dim, tmp2);

  //  printf("B\n");

  mat_mul_mat(dim,      V,  tmp1,      inv_V);
  //  printmat(dim, inv_V);
  mat_mul_mat(dim, 1.0, U0, tmp2, 1.0, inv_V);
  //  printmat(dim, inv_V);
  //  printf("Second GaussJordan\n"); 
  GaussJordan(dim, tmp2, error);
  //  printf("...okay\n");
  mat_mul_mat(dim, inv_V, tmp2, phi);
  //  printmat(dim, phi);

  //  printf("C\n");

#ifndef HAVE_C99
  delete [] inv_V;
  delete [] tmp1;
  delete [] tmp2;
  delete [] T;
  delete [] TN;
  delete [] tmp;
#endif 
}



/* ----------------------------------------------------------------------
 * private method, to get the analytic stiffness coefficients for
 * using the transfer matrix method, zero force boundary condition
 * --------------------------------------------------------------------*/
void greens_function_transfer_matrix_stiffness(int height, int ndof, int dim,
					       double_complex *U0,
					       double_complex *U,
					       double_complex *V,
					       double_complex *phi,
					       Error *error)
{
  //  printf(" \n Entering greens function transfer matrix stiffness \n " ); // dbg_JMM
  int dim_sq = dim*dim;
#ifdef HAVE_C99
  double_complex VT[dim_sq], Vdagger[dim_sq], tmp[dim_sq];
#else
  double_complex *VT, *Vdagger, *tmp;
  VT = new double_complex[dim_sq];
  Vdagger = new double_complex[dim_sq];
  tmp = new double_complex[dim_sq];
#endif

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      /* Note: The factor -1 is the minus from  VT = -V (U + VT)^-1 V^+ */
      MEL(dim, Vdagger, i, j) = -conj(MEL(dim, V, j, i));
    }
  }
  // printf("\n V and Vdagger \n"); // dbg_JMM
  //printmat(dim, V);
  //printmat(dim, Vdagger);

  /* VT = U_N, for which we must flip the order of atoms */
  // reverse_atom_order(ndof, dim, VT, U0); // (ndof, dim, out, in)
  memcpy(VT, U, dim_sq*sizeof(double_complex));

  /* VT is initially -V U^-1 V^+, this corresponds to a fixed boundary */
  GaussJordan(dim, VT, error);
  mat_mul_mat(dim, VT, Vdagger, tmp);
  mat_mul_mat(dim, V, tmp, VT);

  //printf(" \n VT \n"); // dbg_JMM
  //printmat(dim,VT);
  //mat_mul_mat(dim, VT, V, tmp);
  //mat_mul_mat(dim, Vdagger, tmp, VT);

  /* Compute VT */
  if (height < 0) { 
    iterate_Gnn(dim, V, Vdagger, U, VT, 1e-8, 100000, error);
  }
  else if (height > 0) {
    // printf("\n Heightkey iterate Gnn height = %d \n",height); // dbg_JMM
    /* Iterate exactly height */
    iterate_Gnn(dim, V, Vdagger, U, VT, 0.0, height-1, error);   
    //iterate_Gnn(dim, Vdagger, V, U, VT, 0.0, height-1, error);
  }
  if (height != 0) {
    for (int i = 0; i < dim_sq; i++) {
      phi[i] = U0[i] + VT[i];
    }
  } else {
    for (int i = 0; i < dim_sq; i++) {
      phi[i] = U0[i];
    }
  }
  //printf("TAS2 phi \n");printmat(dim, phi);
#ifndef HAVE_C99
  delete [] VT;
  delete [] Vdagger;
  delete [] tmp;
#endif
}



/* ----------------------------------------------------------------------
 * return the displacements within the elastic bulk
 * --------------------------------------------------------------------*/
void transfer_matrix_Gn0(int n, int height, int ndof, int dim,
                         double_complex *U0,
                         double_complex *U,
                         double_complex *V,
                         double_complex *phi,
                         double_complex *Gn0ptr,
                         Error *error)
{
  SquareMatrix<double_complex> Gn0(dim, Gn0ptr);
  SquareMatrix<double_complex> VT(dim), Vdagger(dim);

  /*
   * Compute Green's function and inverse of V
   */
  inverse(Gn0, phi, error);

  /*
   * Compute Vdagger
   */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      /* Note: The factor -1 is the minus from  VT = -V (U + VT)^-1 V^+ */
      Vdagger[i][j] = -conj(MEL(dim, V, j, i));
    }
  }

  for (int i = 0; i < n; i++) {
    SquareMatrix<double_complex> tmp(dim);

    /* VT = U_N, for which we must flip the order of atoms */
    //  reverse_atom_order(ndof, dim, VT, U0);
    VT = U;

    /* VT is initially -V U^-1 V^+, this corresponds to a fixed boundary */
    VT.invert(error);
    mat_mul_mat(VT, Vdagger, tmp);
    mat_mul_mat(V, tmp, VT);

    iterate_Gnn(dim, V, Vdagger.data(), U, VT.data(), 0.0, height-i-1, error);
    VT += U;
    VT.invert(error);

    mat_mul_mat(Vdagger, Gn0, tmp);
    mat_mul_mat(VT, tmp, Gn0);
  }
}



/* ----------------------------------------------------------------------
 * private method, to get the analytic stiffness coefficients for
 * using the renormalization group method, zero force boundary condition
 * --------------------------------------------------------------------*/
void renormalization_group_stiffness(int height, int ndof, int dim,
				     double_complex *U0, double_complex *U,
				     double_complex *V, double_complex *phi,
				     Error *error)
{
  SquareMatrix<double_complex> U1(dim), V1(dim), V1dagger(dim), U01(dim, phi);
  SquareMatrix<double_complex> UN1(dim), inv_U1(dim), tmp1(dim), tmp2(dim);

  int M = lround(log(height)/log(2.0));
  int N = lround(pow(2.0, M));

  if (height < 1 || N != height) {
    char errstr[1024];
    sprintf(errstr, "renormalization_group_stiffness: Height of the system "
	    "needs to be a power of 2 but is %i.", height);
    error->all(FLERR, errstr);
  }

  /*
   * Initialize at first iteration
   */
  U1 = U;
  V1 = V;
  U01 = U0;
  UN1 = U;

  for (int m = 0; m < M; m++) {
    /* V1dagger = V1^+ */
    conj_transpose(V1dagger, V1);
    /* inv_U1 = U1^-1 */
    inverse(inv_U1, U1, error);
    /* tmp1 = inv_U1*V1^+ */
    mat_mul_mat(inv_U1, V1dagger, tmp1);
    /* tmp2 = - V*tmp1, U1 += tmp2, U01 += tmp2 */
    mat_mul_mat(-1.0, V1, tmp1, 0.0, tmp2);
    U1 += tmp2;
    U01 += tmp2;
    /* tmp1 = inv_U1*V1 */
    mat_mul_mat(inv_U1, V1, tmp1);
    /* tmp2 = - V^+*tmp1, U1 += tmp2, UN1 += tmp2 */
    mat_mul_mat(-1.0, V1dagger, tmp1, 0.0, tmp2);
    U1 += tmp2;
    UN1 += tmp2;
    /* V1 = -V1*tmp1 */
    mat_mul_mat(-1.0, V1, tmp1, 0.0, tmp2);
    V1 = tmp2;
  }

  /* V1dagger = V^+ */
  conj_transpose(V1dagger, V1);
  /* inv_U = UN^-1 */
  inverse(inv_U1, UN1, error);
  /* tmp1 = inv_U*V^+ */
  mat_mul_mat(inv_U1, V1dagger, tmp1);
  /* phi = U0 - V*tmp1 */
  U01.mul_then_add(-1.0, V1, tmp1);  
}



/* ----------------------------------------------------------------------
 * symmetrize phi[q=0] such that the total force in x-, y- and z-direction
 * vanishes
 * --------------------------------------------------------------------*/

void enforce_phi0_sum_rule(int dim, double_complex *phi0, double eps)
{
  // assume ndof = 3
  int nlayers = dim/3;

  //  printf("\n enforce nlayers %i \n ",nlayers); // dbg_JMM
  
  double eps2 = eps*eps;

  double sumnorm = 1.0;  // by how much is the sum rule violated?
  double symnorm = 1.0;  // by how much is symmetry violated?
  while (sumnorm > eps2 || symnorm > eps2) {

    // enforce sum rule, i.e. sum of forces in x, y, and z directions
    // should be zero.
    sumnorm = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < dim; j++) {
	double_complex phisum = 0.0;
	for (int k = 0; k < nlayers; k++) {
	  phisum += MEL(dim, phi0, 3*k+i, j);
	  //printf(" phi0(3*k+i,j) %f \n",MEL(dim, phi0, 3*k+i, j)); // dbg_JMM 
	}
	sumnorm += cabs(phisum*conj(phisum));
	//printf("\n sumnorm %f \n",sumnorm); // dbg_JMM
	phisum /= nlayers;
	//printf("\n phisum %f \n",phisum); // dbg_JMM
	for (int k = 0; k < nlayers; k++) {
	  MEL(dim, phi0, j, 3*k+i) -= phisum;
	}
      }
    }

    // symmetrize matrix
    symnorm = 0.0;
    for (int i = 0; i < dim-1; i++) {
      for (int j = i+1; j < dim; j++) {
	double phiavg;
	phiavg = creal(MEL(dim, phi0, i, j)) - creal(MEL(dim, phi0, j, i));
	symnorm += phiavg*phiavg;
	phiavg = 0.5*(creal(MEL(dim, phi0, i, j)) +
		      creal(MEL(dim, phi0, j, i)));
	MEL(dim, phi0, i, j) = phiavg;
	MEL(dim, phi0, j, i) = phiavg;
      }
    }
  }
}

/* --------------------------------------------------------------------*/



} /* namespace */


/*
possible typo at Eqn A15.1 as V may need to oscillate between +/-
mulphase depends on column numbering in geometry file.  each layer needs be set in input file

U + W
- U - W
W + U +
0 W - U

     j-->    j-->
i    U + W +|Y
|    - U - W|- Y
|    W + U +|W + Y
V    - W - U|- W - Y
     Y + W + U + W + Y
 
 
FCC100 is ABABAB. 
Translations (in x dir e.g.) in nn units are 0.0 0.5 0.0 0.5 
(Negative signs possible depending on geometry file indexing.)

FCC111 would be ABCABC
Allowing four layers, for e.g. diamond, 
A B C D A B C D

Shift from lattice plane alph to bet
      bet-->
alph   0  a  b   c   | 0  a   b  c
|     -a  0  b-a c-a |-a  0
|     -b a-b 0   c-b |-b  a-b 0
V     -c a-c b-c 0   |-c

in general, it's shift[j-i]-shift[i] for U loop
and for V loop its -shift[i]+shift[j]

*/
