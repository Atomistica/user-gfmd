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
#ifndef __VEC_H
#define __VEC_H

#include <string.h>

void printvec(int dim, const double_complex *m);
void printvec(int dim, const double *m); 

/*
 * Simple vector class, including support for addition and multiplication
 */

template<typename T, unsigned int stack_dim_ = 0>
class vec {
 public:
    /*!
     * Const copy constructor
     */
    vec(const vec<T, stack_dim_> &other)  : dim_(stack_dim_),
      data_(stack_data_), own_data_(true) {
        if (stack_dim_ == 0) {
            dim_ = other.dim_;
            _alloc();
        }
        assert(dim_ == other.dim_);
        memcpy(data_, other.data_, dim_*sizeof(T));
    }

    /*!
     * Const copy constructor with different dimension
     */
    template<unsigned int other_dim>
    vec(const vec<T, other_dim> &other)  : dim_(stack_dim_),
      own_data_(true), data_(stack_data_) {
        if (stack_dim_ == 0) {
            dim_ = other.dim_;
            _alloc();
        }
        assert(dim_ == other.dim_);
        memcpy(data_, other.data_, dim_*sizeof(T));
    }

#if 0
    /*!
     * Move constructor
     */
    vec(vec<T, stack_dim_> &&other)  : dim_(stack_dim_),
      own_data_(true), data_(stack_data_) {
        if (stack_dim_ == 0) {
            other.own_data_ = false;
            dim_ = other.dim_;
            data_ = other.data;
        }
        else {
            memcpy(data_, other.data_, dim_*sizeof(T));
        }
    }
#endif
 
    vec(const T *data = NULL) : dim_(stack_dim_), data_(stack_data_),
      own_data_(true) {
        assert(stack_dim_ > 0);
        if (data) {
            memcpy(data_, data, dim_*sizeof(T));
        }
    }

    vec(int dim, T *data = NULL, bool copy=false) : dim_(dim),
      data_(stack_data_), own_data_(false) {
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
            std::copy(data, data+dim_, data_);
        }
    }

    ~vec() {
        if (own_data_ && stack_dim_ == 0) {
            delete [] data_;
        }
    }
 
    vec<T, stack_dim_> &operator=(const T &data) {
        fill_with(data);
        return *this;
    }

    /*!
     * Const copy assignment
     */
    vec<T, stack_dim_> &operator=(const vec<T, stack_dim_> &A) {
        set(A.data_);
        return *this;
    }

    /*!
     * Const copy assignment with different dimension
     */
    template<unsigned int Tdim>
    vec<T, stack_dim_> &operator=(const vec<T, Tdim> &A) {
        set(A.data_);
        return *this;
    }

#if 0
    /*!
     * Move assignment
     */
    vec<T, stack_dim_> &operator=(vec<T, stack_dim_> &&A) {
        if (stack_dim_ == 0) {
            data_ = A.data_;
        }
        else {
            set(A.data_);
        }
        return *this;
    }
#endif

    vec<T, stack_dim_> &operator=(const T *data) {
        set(data);
        return *this;
    }

    void set(const T *data) {
        memcpy(data_, data, dim_*sizeof(T));
    }

    T &operator[](int x) {
        return data_[x];
    }

    const T &operator[](int x) const {
        return data_[x];
    }

    vec<T, stack_dim_> &operator+=(const T *A) {
        axpy(1.0, A);
        return *this;
    }

    template<unsigned int Adim>
    vec<T, stack_dim_> &operator+=(const vec<T, Adim> &A) {
        axpy(1.0, A);
        return *this;
    }

    template<typename U>
    vec<T, stack_dim_> operator+(const U &A) const {
        vec<T, stack_dim_> result = *this;
        result += A;
        return result;              
    }

    vec<T, stack_dim_> &operator-=(const T *A) {
        for (int i = 0; i < dim_; i++) {
            data_[i] -= A[i];
        }
        return *this;
    }

    template<unsigned int Adim>
    vec<T, stack_dim_> &operator-=(const vec<T, Adim> &A) {
        for (int i = 0; i < dim_; i++) {
            data_[i] -= A.data_[i];
        }
        return *this;
    }

    template<typename U>
    vec<T, stack_dim_> operator-(const U &A) const {
        vec<T, stack_dim_> result = *this;
        result -= A;
        return result;              
    }

    vec<T, stack_dim_> &operator*=(T alpha) {
        for (int i = 0; i < dim_; i++) {
            data_[i] *= alpha;
        }
        return *this;
    }

    vec<T, stack_dim_> operator*(T alpha) const {
        vec<T, stack_dim_> result = *this;
        result *= alpha;
        return result;              
    }

    vec<T, stack_dim_> &operator/=(T alpha) {
        for (int i = 0; i < dim_; i++) {
            data_[i] /= alpha;
        }
        return *this;
    }

    vec<T, stack_dim_> operator/(T alpha) const {
        vec<T, stack_dim_> result = *this;
        result /= alpha;
        return result;              
    }

    template<typename U, typename V>
    void axpy(U alpha, const V *A) {
        for (int i = 0; i < dim_; i++) {
            data_[i] += alpha*A[i];
        }
    }

    template<typename U, typename V, unsigned int Vdim>
    void axpy(U alpha, const vec<V, Vdim> &A) {
        axpy(alpha, A.const_data());
    }

    void fill_with(T value) {
        for (int i = 0; i < dim_; i++){
            data_[i] = value;
        }
    }

    T *data() {
        return data_;
    }

    const T *const_data() const {
        return data_;
    }

    double nrm2() {
        double acc = 0.0;
        for (int i = 0; i < dim_; i++) {
            acc += data_[i]*data_[i];
        }
        return sqrt(acc);
    }

    double nrm(int ord=2) {
        double acc = 0.0;
        for (int i = 0; i < dim_; i++) {
            acc += pow(data_[i], ord);
        }
        return pow(acc, 1.0/ord);
    }

    bool almost_equal(const T *other, double tol) const {
        bool eq = true;
        for (int i = 0; i < dim_; i++) {
            eq = eq && std::abs(data_[i]-other[i]) < tol;
        }
        return eq;
    }

    template<unsigned int other_dim>
    bool almost_equal(const vec<T, other_dim> &other, double tol) const {
        return almost_equal(other.const_data(), tol);
    }

    void print() {
        printvec(dim_, data_);
    }

    size_t size() const {
        if (stack_dim_ > 0)
            return stack_dim_;
        return dim_;
    }


    int dim_;
    T *data_;
    T stack_data_[stack_dim_];
    bool own_data_;

 protected:
    void _alloc(const T *data = NULL) {
        data_ = new T[dim_];
        own_data_ = true;
        if (data) {
            memcpy(data_, data, dim_*sizeof(T));
        }
        else {
            memset(data_, 0, dim_*sizeof(T));
        }
    }
};


/*!
 * Complex data type need conjugation
 */
template<>
inline double vec<double_complex>::nrm2() {
	double acc = 0.0;
	for (int i = 0; i < dim_; i++) {
		acc += creal(conj(data_[i])*data_[i]);
	}
	return sqrt(acc);
}


/*!
 * Subtract vec from double
 */
template<typename T, unsigned int Tdim>
vec<T, Tdim> operator-(const T *A, const vec<T, Tdim> &B) {
    vec<T, Tdim> result = A;
    result -= B;
    return result;
}


/*!
 * Multiplication with a scalar
 */
template<typename T, unsigned int Tdim>
vec<T, Tdim> operator*(const T &A, const vec<T, Tdim> &B) {
    return B*A;
}


/*!
 * Overload operator<< for stream functions
 */
template<typename T, unsigned int Tdim>
inline std::ostream &operator<<(std::ostream& os, const vec<T, Tdim> &obj)
{
    os << "(";
    for (int i = 0; i < obj.size(); i++) {
        os << obj[i];
        if (i != obj.size()-1)
            os << " ";
    }
    os << ")";

    return os;
}

#endif
