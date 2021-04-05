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
#include <math.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "gtest/gtest.h"

class LinearalgebraTest: public ::testing::Test {
protected:
    LinearalgebraTest() {
    };

    virtual ~LinearalgebraTest() {
    };
    
    /* Tolerance for comparing values */
    static constexpr double tol_ = 1e-12;
};


TEST_F(LinearalgebraTest, dot) {
    double _A[9] = { 1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0 };
    double _x[3] = { 1.0, 2.0, 3.0 };
    double _y[3] = { 4.0, 5.0, 6.0 };

    mat<double> A(3, _A);
    vec<double> x(3, _x);
    vec<double> y(3, _y);

    y = A.dot(x);

    ASSERT_NEAR(y[0], x[0], tol_);
    ASSERT_NEAR(y[1], x[1], tol_);
    ASSERT_NEAR(y[2], x[2], tol_);
}
