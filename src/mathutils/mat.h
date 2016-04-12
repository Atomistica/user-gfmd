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
#ifndef __MAT_H
#define __MAT_H

#include "error.h"
#include "vec.h"

using namespace LAMMPS_NS;

/* backward compatibility */
#define SquareMatrix mat

void printmat(int dim, const double_complex *m);
void printmat(int dim, const double *m); 

#define MEL(dim, m, i, j)  m[(i)*(dim)+(j)]

/*
 * Simple matrix class
 */

template<typename T, unsigned int stack_dim_=0>
class mat {
 public:
    /*!
     * Const copy constructor
     */
    mat(const mat<T, stack_dim_> &other) : dim_(stack_dim_),
      dim_sq_(stack_dim_*stack_dim_), data_(stack_data_), own_data_(true) {
        if (stack_dim_ == 0) {
            dim_ = other.dim_;
            _alloc();
        }
        assert(dim_ == other.dim_);
        memcpy(data_, other.data_, dim_sq_*sizeof(T));
    }

    /*!
     * Const copy constructor with different dimension
     */
    template<unsigned int other_dim>
    mat(const mat<T, other_dim> &other) : dim_(stack_dim_),
      dim_sq_(stack_dim_*stack_dim_), data_(stack_data_), own_data_(true) {
        if (stack_dim_ == 0) {
            dim_ = other.dim_;
            _alloc();
        }
        assert(dim_ == other.dim_);
        memcpy(data_, other.data_, dim_sq_*sizeof(T));
    }

    mat(const T *data = NULL) : dim_(stack_dim_), dim_sq_(stack_dim_*stack_dim_),
      data_(stack_data_), own_data_(true) {
        assert(stack_dim_ > 0);
        if (data) {
            memcpy(data_, data, dim_sq_*sizeof(T));
        }
    }

    mat(int dim, T *data = NULL, bool copy = false) : dim_(dim),
      dim_sq_(dim_*dim_), data_(stack_data_), own_data_(false) {
        if (stack_dim_ == 0) {
            if (data && copy) {
                _alloc(data);
            }
            else if (!data) {
                _alloc();
            }
            else {
                data_ = data;
            }
        }
        else if (data) {
            assert(dim_ == stack_dim_); 
            assert(copy);
            memcpy(data_, data, dim_sq_*sizeof(T));
        }
    }

    mat(int dim, const T *data, bool copy = false) : dim_(dim), dim_sq_(dim_*dim_),
      data_(stack_data_), own_data_(false) {
        assert(data);
        assert(copy);
        if (stack_dim_ == 0) {
            _alloc(data);
        }
        else {
            assert(dim_ == stack_dim_); 
            memcpy(data_, data, dim_sq_*sizeof(T));
        }
    }

    ~mat() {
        if (own_data_ && stack_dim_==0) {
            delete [] data_;
        }
    }
    
    mat<T, stack_dim_> &operator=(const T &data) {
        fill_with(data);    return *this;
    }
    
    /*!
     * Const copy assignment
     */
    mat<T, stack_dim_> &operator=(const mat<T, stack_dim_> &mat) {
        memcpy(data_, mat.data_, dim_sq_*sizeof(T));
        return *this;
    }

    template<typename U>
    mat<T, stack_dim_> &operator=(const U *data) {
        set(data);
        return *this;
    }

    template<typename U>
    void set(const U *data) {
        for (int i = 0; i < dim_sq_; i++) {
            data_[i] = data[i];
        }
    }

    T *operator[](int x) {
        return &MEL(dim_, data_, x, 0);
    }

    const T *operator[](int x) const {
        return &MEL(dim_, data_, x, 0);
    }

    // Matrix addition
    // Return *this for operator chaining 
    // (Guidance from courses.cms.caltech/donnie)
    template<typename U, unsigned int Udim>
    mat<T, stack_dim_> & operator+=(const mat<U, Udim> &A) { 
        axpy(1.0, A);
        return *this;
    }

    template<unsigned int Adim>
    mat<T, stack_dim_> operator+(const mat<T, Adim> &A) const { // no longer const
        mat<T, stack_dim_> result = *this; // New matrix which will be returned
        result += A;                  // Use += to add other to the copy.
        return result;              
    }

    // Marix subtraction
    template<unsigned int Adim>
    mat<T, stack_dim_> & operator-=(const mat<T, Adim> &A) {
        for (int i = 0; i < dim_sq_; i++)
            data_[i] -= A.data_[i];
        return *this;
    }

    template<unsigned int Adim>
    mat<T, stack_dim_> operator-(const mat<T, Adim> &A) const { 
        mat<T> result = *this; // New matrix which will be returned
        result -= A;               // Use -= to subtract other from the copy.
        return result;              
    }

    mat<T, stack_dim_> & operator+=(T *A) {
        for (int i = 0; i < dim_sq_; i++) {
            data_[i] += A[i];
        }
        return *this;
    }
    
    template<typename U, typename V>
    void axpy(U alpha, const V *A) {
        for (int i = 0; i < dim_sq_; i++) {
            data_[i] += alpha*A[i];
        }
    }

    template<typename U, typename V, unsigned int Vdim>
    void axpy(U alpha, const mat<V, Vdim> &A) {
        axpy(alpha, A.const_data());
    }

    template<typename U>
    mat<T, stack_dim_> &operator*=(const U &value) { 
        for (int i = 0; i < dim_sq_; i++) 
            data_[i] *= value;
        return *this;
    }

    mat<T, stack_dim_> operator*(const T &value) const {
        mat<T> result = *this;
        result *= value;
        return result;              
    }

    template<typename U>
    mat<T, stack_dim_> &operator/=(const U &value) { 
        for (int i = 0; i < dim_sq_; i++) 
            data_[i] /= value;
        return *this;
    }

    mat<T, stack_dim_> operator/(const T &value) const {
        mat<T> result = *this;
        result /= value;
        return result;              
    }
    
    mat<T, stack_dim_> &operator-=(T *A) {
        for (int i = 0; i < dim_sq_; i++) {
            data_[i] -= A[i];
        }

        return *this;
    }
    
    vec<T, stack_dim_> dot(const T *A) const {
        vec<T, stack_dim_> result(dim_);
        int m = 0;
        for (int i = 0; i < dim_; i++) {
            result[i] = 0.0;
            for (int j = 0; j < dim_; j++)
                result[i] += data_[m++]*A[j];
        }
        return result;
    }

    vec<T, stack_dim_> dot(const vec<T> &A) const {
        return dot(A.const_data());
    }

    template<typename U>
    void fill_with(U value) {
        for (int i = 0; i < dim_sq_; i++){
            data_[i] = value;
        }
    }

    T *data() {
        return data_;
    }

    const T *const_data() const {
        return data_;
    }

    double norm2() {
        double acc = 0.0;
        for (int i = 0; i < dim_sq_; i++) {
            acc += creal(conj(data_[i])*data_[i]);
        }
        return sqrt(acc);
    }

    double norm(int ord=2) {
        double acc = 0.0;
        for (int i = 0; i < dim_sq_; i++) {
            acc += pow(data_[i], ord);
        }
        return pow(acc, 1.0/ord);
    }
    
    void invert(Error *error) {
        GaussJordan(dim_, data_, error);
    }

    /*
     * Multipy A*B and add to current matrix
     */
    template<unsigned int Adim, unsigned int Bdim>
    void mul_then_add(const mat<T, Adim> &A, const mat<T, Bdim> &B) {
        mat_mul_mat(dim_, 1.0, A.data_, B.data_, 1.0, data_);
    }

    /*
     * Multipy A*B and add to current matrix with prefactors
     */
    template<unsigned int Adim, unsigned int Bdim>
    void mul_then_add(double alpha, const mat<T, Adim> &A,
                      const mat<T, Bdim> &B, double beta=1.0) {
        mat_mul_mat(dim_, alpha, A.data_, B.data_, beta, data_);
    }

    void print() {
        if (stack_dim_ > 0) {
            printmat(stack_dim_, data_);
        }
        else {
            printmat(dim_, data_);
        }
    }

    size_t size() const {
        if (stack_dim_ > 0)
            return stack_dim_;
        return dim_;
    }

    int dim_, dim_sq_;
    T *data_;
    T stack_data_[stack_dim_*stack_dim_];
    bool own_data_;

 protected:
    void _alloc(const T *data = NULL) {
        dim_sq_ = dim_*dim_;
        data_ = new T[dim_sq_];
        own_data_ = true;
        if (data) {
            memcpy(data_, data, dim_sq_*sizeof(T));
        }
        else {
            memset(data_, 0, dim_sq_*sizeof(T));
        }
    }
};


template<typename U, typename T, unsigned int Tdim>
mat<T, Tdim> operator*(const U &A, const mat<T, Tdim> &B) {
    return B*A;
}


/* Overload operator<< for stream functions */
template<typename T, unsigned int Tdim>
inline std::ostream &operator<<(std::ostream& os, const mat<T, Tdim> &obj)
{
    os << "((";
    for (int i = 0; i < obj.size(); i++) {
        for (int j = 0; j < obj.size(); j++) {
            os << obj[i][j];
            if (j != obj.size()-1)
                os << " ";
        }
        if (i != obj.size()-1)
            os << ") (";
        else
            os << ")";
    }
    os << " )";

    return os;
}

#endif
