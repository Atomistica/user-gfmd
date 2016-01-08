#ifndef __LINEARALGEBRA_H
#define __LINEARALGEBRA_H

#include <assert.h>
#include <math.h>

#include "complexcomp.h"
#include "error.h"

using namespace LAMMPS_NS;

// private methods to do matrix operation
// matrix inversion
void GaussJordan(int, double *, Error *error);
void GaussJordan(int, double_complex *, Error *error);

/*!
 * inlined method, to do matrix-matrix multiplication for
 * matric; square matrix is required.
 */
#define MatMulMat mat_mul_mat
template<typename T>
void mat_mul_mat(int dim, double alpha, T *MatA, T *MatB, double beta, T *MatC)
{
  int m=0;
  int idim=0;
  for (int i=0; i<dim; i++){
    for (int j=0; j<dim; j++){
      MatC[m] = beta*MatC[m];
      for (int k=0; k<dim; k++) MatC[m] += alpha*MatA[idim+k]*MatB[k*dim+j];
      m++;
    }
    idim += dim;
  }
}


/*!
 * inlined method, to do matrix-matrix multiplication for
 * matric; square matrix is required.
 */
#define MatMulMat mat_mul_mat
template<typename T>
void mat_mul_mat(int dim, T *MatA, T *MatB, T *MatC)
{
  mat_mul_mat(dim, 1.0, MatA, MatB, 0.0, MatC);
}

/* vector and matrix helper classes */
#include "vec.h"
#include "mat.h"

// assemble small matrices into large one
void put_matrix(int, double_complex *, int, double_complex *, int, int);
void put_conj_matrix(int, double_complex *, int, double_complex *, int, int);
void take_matrix(int, double_complex *, int, double_complex *, int, int);

// to get the sorted eigenvalue/eigenvectors of a real,symmetric matrix
int  _jacobi(double **, int, double *, double **);
void _rotate(double **, int, int, int, int, double, double);

bool _is_Hermitian(int, double_complex *);
bool _is_Hermitian(int, double_complex *, double tol);

#ifdef GFMD_CHECK_POSITIVE_DEFINITE
void eigvals(int, double *, double *, Error *);
void eigvals(int, double_complex *, double *, Error *);

template<typename T>
bool is_positive_definite(int n, T *m, double *w, Error *error)
{
    double _w[n];
    if (!w)  w = _w;
    eigvals(n, m, w, error);
    for (int k = 0; k < n; k++) {
        if (w[k] <= 0.0)
            return false;
    }

    return true;
}

template<typename T>
void warn_positive_definite(char *name, int n, T *m, Error *error)
{
    double w[n];
    if (!is_positive_definite(n, m, w, error)) {
        char errstr[1024];
        sprintf(errstr, "Matrix %s is not positive definite. "
                "Eigenvalues are %f", name, w[0]);
        for (int l = 1; l < n; l++) {
            sprintf(errstr, "%s, %f", errstr, w[l]);
        }
        //error->all(FLERR,errstr);
        printf("%s\n", errstr);
        printf("Offending matrix:\n", name);
        printmat(n, m);
    }
}

#endif

void unit_outer_complex(int, double *, double *, double_complex *); 
double vec_norm2(int, double *);


/* ----------------------------------------------------------------------
 * inlined method, invert a 3x3 matrix
 * --------------------------------------------------------------------*/
template<typename T>
void invert3x3(T *A)
{
  T Ainv[9];

#define I3(i, j)  ((j-1)*3+(i-1))

  Ainv[I3(1,1)] = A[I3(2,2)]*A[I3(3,3)] - A[I3(3,2)]*A[I3(2,3)];
  Ainv[I3(1,2)] = A[I3(3,2)]*A[I3(1,3)] - A[I3(1,2)]*A[I3(3,3)];
  Ainv[I3(1,3)] = A[I3(1,2)]*A[I3(2,3)] - A[I3(1,3)]*A[I3(2,2)];
  T detA = Ainv[I3(1,1)]*A[I3(1,1)] + Ainv[I3(1,2)]*A[I3(2,1)] +
    Ainv[I3(1,3)]*A[I3(3,1)];

  Ainv[I3(1,1)] = Ainv[I3(1,1)]/detA;
  Ainv[I3(1,2)] = Ainv[I3(1,2)]/detA;
  Ainv[I3(1,3)] = Ainv[I3(1,3)]/detA;

  Ainv[I3(2,1)] = (A[I3(2,3)]*A[I3(3,1)] - A[I3(2,1)]*A[I3(3,3)])/detA;
  Ainv[I3(2,2)] = (A[I3(1,1)]*A[I3(3,3)] - A[I3(3,1)]*A[I3(1,3)])/detA;
  Ainv[I3(2,3)] = (A[I3(2,1)]*A[I3(1,3)] - A[I3(1,1)]*A[I3(2,3)])/detA;

  Ainv[I3(3,1)] = (A[I3(2,1)]*A[I3(3,2)] - A[I3(2,2)]*A[I3(3,1)])/detA;
  Ainv[I3(3,2)] = (A[I3(3,1)]*A[I3(1,2)] - A[I3(1,1)]*A[I3(3,2)])/detA;
  Ainv[I3(3,3)] = (A[I3(1,1)]*A[I3(2,2)] - A[I3(1,2)]*A[I3(2,1)])/detA;

  for (int i = 0; i < 9; i++)
    A[i] = Ainv[i];

#undef I3
}

/* ----------------------------------------------------------------------
 * inlined method, to do vector dot product
 * --------------------------------------------------------------------*/
template<typename T>
T vec_dot_vec(int dim, T *v1, T *v2)
{
  T dot = 0.0;
  for (int i = 0; i < dim; i++) {
    dot += v1[i]*v2[i];
  }
  return dot;
}


/* ----------------------------------------------------------------------
 * inlined method, to do vector dot product
 * --------------------------------------------------------------------*/
template<typename T>
T conj_vec_dot_vec(int dim, T *v1, T *v2)
{
  T dot = 0.0;
  for (int i = 0; i < dim; i++) {
    dot += conj(v1[i])*v2[i];
  }
  return dot;
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for
 * matrix and vector; square matrix is required and output vector is
 * overwritten.
 * --------------------------------------------------------------------*/
#define MatMulVec mat_mul_vec
template<typename T>
void mat_mul_vec(int dim, const T *Mat, const T *Vin, T *Vout)
{
  for (int i = 0; i < dim; i++, Vout++){
    *Vout = 0.0;
    for (int j = 0; j < dim; j++, Mat++) {
      *Vout += (*Mat)*Vin[j];
    }
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for conjugate of a
 * cmatrix and vector; square matrix is required and output vector is
 * overwritten.
 * --------------------------------------------------------------------*/
template<typename T>
void conj_mat_mul_vec(int dim, T *Mat, T *Vin, T *Vout)
{
  for (int i = 0; i < dim; i++){
    int m = i;
    Vout[i] = 0.0;
    for (int j = 0; j < dim; j++, m+=dim)
      Vout[i] += conj(Mat[m])*Vin[j];
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for
 * matrix and vector; square matrix is required and result is added to
 * the output vector.
 * --------------------------------------------------------------------*/
#define MatMulAddVec mat_mul_add_vec
template<typename T>
void mat_mul_add_vec(int dim, T *Mat, T *Vin, T *Vout)
{
  int m=0;
  for (int i=0; i<dim; i++){
    for (int j=0; j<dim; j++) Vout[i] += Mat[m++]*Vin[j];
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for matrix and
 * vector; square matrix is required and result is subtracted from the
 * output vector.
 * --------------------------------------------------------------------*/
#define MatMulSubVec mat_mul_sub_vec
template<typename T>
void mat_mul_sub_vec(int dim, T *Mat, T *Vin, T *Vout)
{
  int m=0;
  for (int i=0; i<dim; i++){
    for (int j=0; j<dim; j++) Vout[i] -= Mat[m++]*Vin[j];
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for the conjugate
 * matrix and vector; square matrix is required and result is added to
 * the output vector.
 * --------------------------------------------------------------------*/
template<typename T>
void conj_mat_mul_add_vec(int dim, T *Mat, T *Vin, T *Vout)
{
  for (int i = 0; i < dim; i++){
    int m = i;
    for (int j = 0; j < dim; j++, m+=dim)
      Vout[i] += conj(Mat[m])*Vin[j];
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-vector multiplication for the conjugate
 * matrix and vector; square matrix is required and result is subtracted
 * from the output vector.
 * --------------------------------------------------------------------*/
template<typename T>
void conj_mat_mul_sub_vec(int dim, T *Mat, T *Vin, T *Vout)
{
  for (int i = 0; i < dim; i++) {
    int m = i;
    for (int j = 0; j < dim; j++, m+=dim)
      Vout[i] -= conj(Mat[m])*Vin[j];
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-matrix multiplication for
 * matric; square matrix is required.
 * --------------------------------------------------------------------*/
template<typename T>
void conj_mat_mul_mat(int dim, double alpha, T *MatA, T *MatB, double beta,
		      T *MatC)
{
  int m=0;
  for (int i=0; i<dim; i++){
    for (int j=0; j<dim; j++){
      MatC[m] = beta*MatC[m];
      for (int k=0; k<dim; k++)
	MatC[m] += alpha*conj(MatA[k*dim+i])*MatB[k*dim+j];
      m++;
    }
  }
}


/* ----------------------------------------------------------------------
 * inlined method, to do matrix-matrix multiplication for
 * matric; square matrix is required.
 * --------------------------------------------------------------------*/
template<typename T>
void conj_mat_mul_mat(int dim, T *MatA, T *MatB, T *MatC)
{
  conj_mat_mul_mat(dim, 1.0, MatA, MatB, 0.0, MatC);
}


/*
 * Multiply A*B and override output matrix
 */
template<typename T>
void mat_mul_mat(const SquareMatrix<T> &MatA, const SquareMatrix<T> &MatB,
		 SquareMatrix<T> &MatC)
{
  mat_mul_mat(MatA.dim_, 1.0, MatA.data_, MatB.data_, 0.0, MatC.data_);
}


/*
 * Multiply A*B and override output matrix
 */
template<typename T>
void mat_mul_mat(double_complex *MatA, const SquareMatrix<T> &MatB,
                 SquareMatrix<T> &MatC)
{
  mat_mul_mat(MatB.dim_, 1.0, MatA, MatB.data_, 0.0, MatC.data_);
}


/*
 * Multiply A*B and override output matrix
 */
template<typename T>
void mat_mul_mat(const SquareMatrix<T> &MatA, double_complex *MatB,
                 SquareMatrix<T> &MatC)
{
  mat_mul_mat(MatA.dim_, 1.0, MatA.data_, MatB, 0.0, MatC.data_);
}


/*
 * Multiply A*B and override output matrix
 */
template<typename T>
void mat_mul_mat(double_complex *MatA, double_complex *MatB,
                 SquareMatrix<T> &MatC)
{
  mat_mul_mat(MatC.dim_, 1.0, MatA, MatB, 0.0, MatC.data_);
}


/*
 * Multiply A*B and override output matrix
 */
template<typename T>
void mat_mul_mat(double alpha, const SquareMatrix<T> &MatA,
		 const SquareMatrix<T> &MatB, double beta,
		 SquareMatrix<T> &MatC)
{
  mat_mul_mat(MatA.dim_, alpha, MatA.data_, MatB.data_, beta, MatC.data_);
}


/*!
 * Hermitian adjoint
 */
template<typename T>
inline void conj_transpose(int dim, T *mat1, const T *mat2)
{
  int m = 0;
  for (int i = 0; i < dim; i++){
    int n = i;
    for (int j = 0; j < dim; j++, m++, n += dim)
      mat1[m] = mat2[n];
  }
}


/*!
 * Hermitian adjoint, complex type need conjugation
 */
template<>
inline void conj_transpose<double_complex>(int dim, double_complex *mat1,
                                           const double_complex *mat2)
{
  int m = 0;
  for (int i = 0; i < dim; i++){
    int n = i;
    for (int j = 0; j < dim; j++, m++, n += dim)
      mat1[m] = conj(mat2[n]);
  }
}


/*
 * Hermitian adjoint
 */
template<typename T>
inline void conj_transpose(mat<T> &res, const mat<T> &val)
{
  conj_transpose(res.dim_, res.data_, val.data_);
}


/*!
 * Compute the inverse
 */
template<typename T>
void inverse(SquareMatrix<T> &res, const SquareMatrix<T> &val, Error *error)
{
  res = val;
  res.invert(error);
}


/*!
 * Compute the inverse
 */
template<typename T>
void inverse(SquareMatrix<T> &res, double *val, Error *error)
{
  res = val;
  res.invert(error);
}


/*!
 * Compute the inverse
 */
template<typename T>
void inverse(SquareMatrix<T> &res, double_complex *val, Error *error)
{
  res = val;
  res.invert(error);
}

/*!
 * Identity matrix
 */
template<typename T>
mat<T> identity(int n) {
	mat<T> id(n);
	for (int i = 0; i < n; i++)
		id[i][i] = 1;
	return id;
}

/*!
 * Outer product
 */
template<typename T, unsigned int Adim, unsigned int Bdim, unsigned int Cdim>
void outer(const vec<T, Adim> &A, const vec<T, Bdim> &B, mat<T, Cdim> &C)
{
	assert(A.dim_ == B.dim_);
	assert(A.dim_ == C.dim_);

	for (int i = 0; i < A.dim_; i++) {
		for (int j = 0; j < A.dim_; j++) {
			C[i][j] = A[i]*B[j];
		}
	}
}


/*!
 * Outer product
 */
template<typename T, unsigned int Adim, unsigned int Bdim>
mat<T> outer(const vec<T, Adim> &A, const vec<T, Bdim> &B)
{
    mat<T> C(A.dim_);
    outer(A, B, C);

    return C;
}


double string_to_double(const char *file, int line, Error *error, char *str);
void read_double(Error *error, int narg, char **arg, int *carg, double *ptr,
                 const char *description, double rangemin, double rangemax,
                 int prnt=0);
void read_integer(Error *error, int narg, char **arg, int *carg, int *ptr,
                  const char *description, int rangemin, int rangemax,
                  bool prnt=false);

inline int modulo(int i, int n)
{
    while (i < 0) i += n;
    return i%n;
}

#endif
