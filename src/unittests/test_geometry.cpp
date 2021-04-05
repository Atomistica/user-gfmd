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

#include "geometry.h"

#include "gtest/gtest.h"

class GeometryTest: public ::testing::Test {
protected:
    GeometryTest() {
    };

    virtual ~GeometryTest() {
    };
    
    /* Tolerance for comparing values */
    static constexpr double tol_ = 1e-12;
};


TEST_F(GeometryTest, clockwise_angle) {
    double alpha = 180*clockwise_angle(1.0,0.0, 0.0,1.0)/M_PI;
    ASSERT_NEAR(alpha, -90.0, tol_);
    
    alpha = 180*clockwise_angle(0.0,1.0, 1.0,0.0)/M_PI;
    ASSERT_NEAR(alpha, 90.0, tol_);
    
    alpha = 180*clockwise_angle(1.0,0.0, -1.0,1.0)/M_PI;
    ASSERT_NEAR(alpha, -135.0, tol_);
}


TEST_F(GeometryTest, lines_cross) {
    ASSERT_TRUE(lines_cross(0,0, 1,0,  0.5,-0.5, 0.5,0.5));
    ASSERT_FALSE(lines_cross(0,0, 1,0,  0.5,0.5, 0.5,1.5));
    
    ASSERT_TRUE(lines_cross(0.5,-0.5, 0.5,0.5,  0,0, 1,0));
    ASSERT_FALSE(lines_cross(0.5,0.5, 0.5,1.5,  0,0, 1,0));
}


TEST_F(GeometryTest, tr_is_inside) {
    ASSERT_TRUE(tr_is_inside(0,0, -2,1, 1,-2, 1,1));
    ASSERT_TRUE(tr_is_inside(0,0, 1,-2, -2,1, 1,1));
    
    ASSERT_FALSE(tr_is_inside(2,2, 1,-2, -2,1, 1,1));
    ASSERT_FALSE(tr_is_inside(-3,0, 1,-2, -2,1, 1,1));
}


TEST_F(GeometryTest, tr_area) {
    ASSERT_NEAR(tr_area(-2,1, 1,-2, 1,1), 4.5, tol_);
}


TEST_F(GeometryTest, tr_overlap) {
    ASSERT_TRUE(tr_overlap(-2,1, -2,0, 1,-2, 1,1));
    ASSERT_FALSE(tr_overlap(-2,1, 2,0, 1,-2, 1,1));
}



