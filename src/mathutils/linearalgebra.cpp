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
#include <string.h>
#include <stdlib.h>

#include "linearalgebra.h"

#ifdef HAVE_MKL
#include "mkl_lapack.h"
#endif

#define MAXJACOBI 100

void put_matrix(int ldim, double_complex *lmat, int sdim, double_complex *smat,
                int y, int x)
{
  for (int j = 0; j < sdim; j++) {
    for (int i = 0; i < sdim; i++) {
      MEL(ldim, lmat, y+j, x+i) = MEL(sdim, smat, j, i);
    }
  }
}


void put_conj_matrix(int ldim, double_complex *lmat, int sdim,
                     double_complex *smat, int y, int x)
{
  for (int j = 0; j < sdim; j++) {
    for (int i = 0; i < sdim; i++) {
      MEL(ldim, lmat, y+j, x+i) = conj(MEL(sdim, smat, i, j));
    }
  }
}


void take_matrix(int ldim, double_complex *lmat, int sdim,
                 double_complex *smat, int y, int x)
{
  for (int j = 0; j < sdim; j++) {
    for (int i = 0; i < sdim; i++) {
      MEL(sdim, smat, j, i) = MEL(ldim, lmat, y+j, x+i);
    }
  }
}


bool _is_Hermitian(int dim, double_complex *mat)
{
  bool Hermitian = true;

  int m = 0;
  for (int i = 0; i < dim; i++){
    int n = i;
    for (int j = 0; j < dim; j++, m++, n+=dim)
      Hermitian = Hermitian && ( cabs(mat[m] - conj(mat[n])) < 1e-6 );
  }

  return Hermitian;
}


bool _is_Hermitian(int dim, double_complex *mat, double tol)
{
  bool Hermitian = true;
  double_complex trace = 0.0;
  for (int i = 0; i < dim; i++) trace = trace + mat[i*dim+i];

  int m = 0;
  for (int i = 0; i < dim; i++){
    int n = i;
    for (int j = 0; j < dim; j++, m++, n+=dim) {
      
      if (Hermitian && ( cabs(mat[m] - conj(mat[n]))/cabs(trace) > tol )) { 
        printf("%f %f  != %f %f \n", creal(mat[m]),cimag(mat[m]),creal(mat[n]),cimag(mat[n]));
	printf("at %d %d and %d %d ( %d mod %d )\n", m/dim, m % dim, n/dim, n % dim, n, dim);
        printf("  because  | %f, %f | / %f  = %f > tolerance (%f) ", creal(mat[m] - conj(mat[n])),cimag(mat[m] - conj(mat[n])),cabs(trace),cabs(mat[m] - conj(mat[n]))/cabs(trace), tol);
      }

      Hermitian = Hermitian && ( cabs(mat[m] - conj(mat[n]))/cabs(trace) <= tol );
    }
  }

  return Hermitian;
}



/* ----------------------------------------------------------------------
 * Private method, to compute all eigenvalues and eigenvectors of a real
 * symmetric matrix a[0..n-1][0..n-1]. On output, elements of a above
 * the diagonal are destroyed. d[0..n-1] returns the eigenvalues of a.
 * v[0..n-1][0..n-1] is a matrix whose columns contain, on output, the
 * normalized eigenvectors of a.
 * If maximum iteration is reached it exits and returns 1 instead of 0.
 *
 * Adapted from the Numerical Recipes in Fortran
 * --------------------------------------------------------------------*/

int _jacobi(double **matrix, int n, double *evalues, double **evectors)
{
  int i,j,k;
  double tresh,theta,tau,t,sm,s,h,g,c;
  double *b, *z;

  b = new double [n];
  z = new double [n];

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) evectors[i][j] = 0.;
    evectors[i][i] = 1.0;
  }
  for (i = 0; i < n; i++) {
    b[i] = evalues[i] = matrix[i][i];
    z[i] = 0.;
  }

  for (int iter = 1; iter <= MAXJACOBI; iter++) {
    sm = 0.0;
    for (i = 0; i < n-1; i++)
      for (j = i+1; j < n; j++) sm += fabs(matrix[i][j]);
    if (sm == 0.0){
      for (i=0;i<n-1;i++) { // sort eigen value and eigen vectors
        double p = evalues[k=i];
        for (j=i+1;j<n;j++) if (evalues[j] < p) p = evalues[k=j];
        if (k != i) {
          evalues[k] = evalues[i];
          evalues[i] = p;
          for (j=0;j<n;j++) {
            p = evectors[j][i];
            evectors[j][i] = evectors[j][k];
            evectors[j][k] = p;
          }
        }
      }

      delete []b;
      delete []z;

      return 0;
    }

    if (iter < 4) tresh = 0.2*sm/(n*n);
    else tresh = 0.0;

    for (i = 0; i < n-1; i++) {
      for (j = i+1; j < n; j++) {
        g = 100.0*fabs(matrix[i][j]);
        if (iter > 4 && fabs(evalues[i])+g == fabs(evalues[i])
        && fabs(evalues[j])+g == fabs(evalues[j]))
          matrix[i][j] = 0.0;
        else if (fabs(matrix[i][j]) > tresh) {
          h = evalues[j]-evalues[i];
          if (fabs(h)+g == fabs(h)) t = (matrix[i][j])/h;
          else {
            theta = 0.5*h/(matrix[i][j]);
            t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c = 1.0/sqrt(1.0+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*matrix[i][j];
          z[i] -= h;
          z[j] += h;
          evalues[i] -= h;
          evalues[j] += h;
          matrix[i][j] = 0.0;
          for (k = 0;   k < i; k++) _rotate(matrix,k,i,k,j,s,tau);
          for (k = i+1; k < j; k++) _rotate(matrix,i,k,k,j,s,tau);
          for (k = j+1; k < n; k++) _rotate(matrix,i,k,j,k,s,tau);
          for (k = 0;   k < n; k++) _rotate(evectors,k,i,k,j,s,tau);
        }
      }
    }

    for (i = 0; i < n; i++) {
      evalues[i] = b[i] += z[i];
      z[i] = 0.0;
    }
  }
  delete []b;
  delete []z;
  return 1;
}

/* ----------------------------------------------------------------------
 * perform a single Jacobi rotation
 * ------------------------------------------------------------------- */

void _rotate(double **matrix, int i, int j, int k, int l, double s, double tau)
{
  double g = matrix[i][j];
  double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}


#ifdef HAVE_MKL

void GaussJordan(int n, double *Mat, Error *error)
{
  int n_sq = n*n;

  int ipiv[n];
  int info;

  double tmp1[n_sq], tmp2[n_sq];

  memcpy(tmp1, Mat, n_sq*sizeof(double));
  memset(tmp2, 0, n_sq*sizeof(double));
  for (int i = 0; i < n; i++)
    MEL(n, tmp2, i, i) = 1.0;

  dgesv(&n, &n, tmp1, &n, ipiv, tmp2, &n, &info);

  if (info)
    error->all("dgesv failed.");

  memcpy(Mat, tmp2, n_sq*sizeof(double));
}


void GaussJordan(int n, double_complex *Mat, Error *error)
{
  int n_sq = n*n;

  int ipiv[n];
  int info;

  double_complex tmp1[n_sq], tmp2[n_sq];

  memcpy(tmp1, Mat, n_sq*sizeof(double_complex));
  memset(tmp2, 0, n_sq*sizeof(double_complex));
  for (int i = 0; i < n; i++)
    MEL(n, tmp2, i, i) = 1.0;

  zgesv(&n, &n, (MKL_Complex16*) tmp1, &n, ipiv, (MKL_Complex16*) tmp2,
        &n, &info);

  if (info)
    error->all("zgesv failed.");

  memcpy(Mat, tmp2, n_sq*sizeof(double_complex));
}

#else

/* ----------------------------------------------------------------------
 * private method, to get the inverse of a double precision matrix
 * by means of Gaussian-Jordan Elimination with full pivoting.
 *
 * Adapted from the Numerical Recipes in Fortran.
 * --------------------------------------------------------------------*/
void GaussJordan(int n, double *Mat, Error *error)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
#ifdef HAVE_C99
  int indxc[n],indxr[n],ipiv[n];
#else
  int *indxc, *indxr, *ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
#endif
  double big, dum, pivinv;

  for (i=0; i<n; i++) ipiv[i] = 0;
  for (i=0; i<n; i++){
    big = 0.;
    for (j=0; j<n; j++){
      if (ipiv[j] != 1){
        for (k=0; k<n; k++){
          if (ipiv[k] == 0){
            idr = j*n+k;
            if (fabs(Mat[idr]) >= big){
              big  = fabs(Mat[idr]);
              irow = j;
              icol = k;
            }
          }else if (ipiv[k] >1){ 
            error->one(FLERR,"FixGFMD: Singular matrix in double GaussJordan!");
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l=0; l<n; l++){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == 0.) error->one(FLERR,"FixGFMD: Singular matrix in double GaussJordan!");
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.;
    idr = icol*n;
    for (l=0; l<n; l++) Mat[idr+l] *= pivinv;
    for (ll=0; ll<n; ll++){
      if (ll != icol){
        idc = ll*n+icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l=0; l<n; l++) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k=0; k<n; k++){
        idr = k*n+rl;
        idc = k*n+cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }
#ifndef HAVE_C99
  delete [] indxr;
  delete [] indxc;
  delete [] ipiv;
#endif
return;
}

/* ----------------------------------------------------------------------
 * private method, to get the inverse of a complex matrix by means of
 * Gaussian-Jordan Elimination with full pivoting.
 *
 * Adapted from the Numerical Recipes in Fortran.
 * --------------------------------------------------------------------*/
void GaussJordan(int n, double_complex *Mat, Error *error)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
#ifdef HAVE_C99
  int indxc[n],indxr[n],ipiv[n];
#else
  int *indxc, *indxr, *ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
#endif
  double big, nmjk;
  double_complex dum, pivinv;

  for (i=0; i<n; i++) ipiv[i] = 0;
  for (i=0; i<n; i++){
    big = 0.;
    for (j=0; j<n; j++){
      if (ipiv[j] != 1){
        for (k=0; k<n; k++){
          if (ipiv[k] == 0){
            idr = j*n+k;
            nmjk = cabs(Mat[idr]);
            if (nmjk >= big){
              big  = nmjk;
              irow = j;
              icol = k;
            }
          }else if (ipiv[k]>1){ printmat(n, Mat);
            error->one(FLERR,"FixGFMD: Singular matrix in complex GaussJordan!");
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l=0; l<n; l++){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == 0.0) { printmat(n, Mat);
      error->one(FLERR,"FixGFMD: Singular matrix in complex GaussJordan!");}

    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.0;
    idr = icol*n;
    for (l=0; l<n; l++) Mat[idr+l] *= pivinv;
    for (ll=0; ll<n; ll++){
      if (ll != icol){
        idc = ll*n+icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l=0; l<n; l++) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k=0; k<n; k++){
        idr = k*n+rl;
        idc = k*n+cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }
#ifndef HAVE_C99
  delete [] indxr;
  delete [] indxc;
  delete [] ipiv;
#endif
return;
}

#endif  

/* ----------------------------------------------------------------------
 * Outer product of two normalized vectors
 * --------------------------------------------------------------------*/
// For complex vectors, also conjugate v2 conventionally
// memory for outer must alrady by allocated
void unit_outer_complex(int dim, double *v1, double *v2, double_complex *outer)
{

  double magsqv1 = 0.0, magv1;
  double magsqv2 = 0.0, magv2;

  memset(outer, 0, dim*dim*sizeof(double_complex));

  for (int i = 0; i < dim; i++) {
    magsqv1 += pow(v1[i], 2.0);
    magsqv2 += pow(v2[i], 2.0);
  }
  magv1 = sqrt(magsqv1);
  magv2 = sqrt(magsqv2);

  // TAS i*dim + j is row major order
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      outer[i*dim + j] = (double_complex)(v1[i]/magv1*v2[j]/magv2);
}

/* ----------------------------------------------------------------------
 * Vector norm
 * --------------------------------------------------------------------*/

double vec_norm2(int dim, double *v1)
{
  double norm2 = 0.0;
  for (int i = 0; i < dim; i++) 
    norm2 += (v1[i] * v1[i]);
  return norm2;
}






/* ----------------------------------------------------------------------
 * print a complex vector to screen
 * --------------------------------------------------------------------*/

void printvec(int dim, const double_complex *m)
{

  char str[1024];

  bool real_only = false;

  for (int i = 0; i < dim; i++) {
    if (fabs(cimag(m[i])) > 1e-12)
      real_only = false;
  }

  strcpy(str, " { ");

  for (int l = 0; l < dim; l++) {
    if (real_only) {
      sprintf(str, "%s%10.3e", str, creal(m[l]));
    }
    else {
      sprintf(str, "%s%10.3e+i*%10.3e", str, creal(m[l]), cimag(m[l]));
    }
    if (l != dim-1)  sprintf(str, "%s, ", str);
  }

  printf("%s },\n", str);
}

/* ----------------------------------------------------------------------
 * print a real vector/ matrix to screen
 * --------------------------------------------------------------------*/
#if 1 
void printvec(int dim, const double *m)
{
  char str[1024];
  strcpy(str, " { ");
  for (int l = 0; l < dim; l++) {
    sprintf(str, "%s%10.3e", str, m[l]);
    if (l != dim-1)  sprintf(str, "%s, ", str);
  }
  printf("%s },\n", str);
}
void printmat(int dim, const double *m)
{
  int dim_sq = dim*dim;
  char str[1024];

  for (int k = 0; k < dim; k++) {
    if (k == 0) 
      strcpy(str, "{{ ");
    else
      strcpy(str, " { ");

    for (int l = 0; l < dim; l++) {
      //printf("%d %d",k,l);
      sprintf(str, "%s%10.3e", str, (float)(MEL(dim, m, k, l)));
      if (l != dim-1)  sprintf(str, "%s, ", str);
    }
    if (k == dim-1)
      printf("%s }}\n", str);
    else
      printf("%s },\n", str);
  }
}

#endif

/* ----------------------------------------------------------------------
 * print a complex matrix to screen
 * --------------------------------------------------------------------*/

void printmat(int dim, const double_complex *m)
{
  int dim_sq = dim*dim;
  char str[1024];

  bool real_only = false;

  for (int i = 0; i < dim_sq; i++) {
    if (fabs(cimag(m[i])) > 1e-12)
      real_only = false;
  }

  for (int k = 0; k < dim; k++) {
    if (k == 0) 
      strcpy(str, "{{ ");
    else
      strcpy(str, " { ");

    for (int l = 0; l < dim; l++) {
      if (real_only) {
	sprintf(str, "%s%10.3e", str, creal(MEL(dim, m, k, l)));
      }
      else {
	  sprintf(str, "%s%10.3e+i*%10.3e", str, creal(MEL(dim, m, k, l)),
		  cimag(MEL(dim, m, k, l)));
      }
      if (l != dim-1)  sprintf(str, "%s, ", str);
    }

    if (k == dim-1)
      printf("%s }}\n", str);
    else
      printf("%s },\n", str);
  }
}

double string_to_double(const char *file, int line, Error *error, char *str)
{
  int n = strlen(str);
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    error->all(file,line,"Expected floating point parameter "
               "in input script or data file");
  }

  return atof(str);
}

void read_double(Error *error, int narg, char **arg, int *carg, double *ptr,
                 const char *description, double rangemin, double rangemax,
                 int prnt) {
    if (prnt) printf(description);

    if (*carg >= narg)
        error->all(FLERR,"Missing input value.\n");
    
    if (prnt) puts(arg[*carg]);
    char *endptr;
    *ptr = strtod(arg[*carg], &endptr);
    
    char errstr[1024];
    if (endptr == arg[*carg]) {
        sprintf(errstr, "Error converting to number.");
        error->all(FLERR, errstr);
    }
    if (rangemin != rangemax) {
        if ((*ptr < rangemin) || (*ptr > rangemax)) {
            sprintf(errstr, "Must be greater than %e and less than %e", 
                    rangemin, rangemax);
            error->all(FLERR, errstr);
        }
    }

    (*carg)++;
}


// This is a near-identical copy of read_dbl and could be combined with templates
void read_integer(Error *error, int narg, char **arg, int *carg, int *ptr,
                  const char *description, int rangemin, int rangemax,
                  bool prnt)
{
    if (prnt) printf(description);
    
    if (*carg >= narg)
        error->all(FLERR,"Missing input value.\n");

    if (prnt) puts(arg[*carg]);
    char *endptr;
    *ptr = strtol(arg[*carg], &endptr, 10);
    
    char errstr[1024];
    if (endptr == arg[*carg]) {
        sprintf(errstr, "Error converting to number.");
        error->all(FLERR, errstr);
    }
    if (rangemin != rangemax) {
        if ((*ptr < rangemin) || (*ptr > rangemax)) {
            sprintf(errstr, "Must be greater than %i and less than %i",
                    rangemin, rangemax);
            error->all(FLERR, errstr);
        }
    }
    
    (*carg)++;
}

#if GFMD_CHECK_POSITIVE_DEFINITE
extern "C"
void dsyev_(const char *, const char *, int *, double *, int *, double *, double *, int *, int *);
extern "C"
void zheev_(const char *, const char *, int *, double_complex *, int *, double *, double_complex *, int *, double *, int *);

void eigvals(int n, double_complex *m, double *w, Error *error)
{
    double_complex A[n*n];
    int lwork = std::max(1, 2*n-1);
    double_complex work[lwork];
    int lrwork = std::max(1, 3*n-2);
    double rwork[lrwork];
    std::copy(m, m+n*n, A);
    int info;
    zheev_("N", "U", &n, A, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
        error->one(FLERR,"zheev failed.");
    }
}

void eigvals(int n, double *m, double *w, Error *error)
{
    double A[n*n];
    int lwork = std::max(1, 3*n-1);
    double work[lwork];
    std::copy(m, m+n*n, A);
    int info;
    dsyev_("N", "U", &n, A, &n, w, work, &lwork, &info);
    if (info != 0) {
        error->one(FLERR,"dsyev failed.");
    }
}

#endif
