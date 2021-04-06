/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

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

#include "boussinesq_cerruti.h"
#include "geometry.h"
#include "li_berger.h"

#include "gtest/gtest.h"

using namespace LiBerger;

class LiBergerTest: public ::testing::Test {
protected:
    LiBergerTest() {
    };

    virtual ~LiBergerTest() {
    };
    
    static constexpr int n_trials_ = 100;
    static constexpr double G_ = 1.0;
    static constexpr double nu_ = 1.0/3;

    /* Tolerance for comparing values */
    static constexpr double tol_ = 1e-12;
};


TEST_F(LiBergerTest, check_triangle_solution_symmetry) {
    /*
     * Check the symmetry of the Green's function
     */
    for (int i = 0; i < n_trials_; i++) {
        double A1x = 2*rand1()-1, A1y = 2*rand1()-1;
        double A2x = 2*rand1()-1, A2y = 2*rand1()-1;
        double A3x = 2*rand1()-1, A3y = 2*rand1()-1;
        double Bx = 10*(2*rand1()-1), By = 10*(2*rand1()-1);
       
        double uxx = 0.0, uyx = 0.0, uzx = 0.0;
        tr_uniform_PQ(G_,nu_, 0,1,0, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, uxx,uyx,uzx);

        double uxy = 0.0, uyy = 0.0, uzy = 0.0;
        tr_uniform_PQ(G_,nu_, 0,0,1, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, uxy,uyy,uzy);

        double uxz = 0.0, uyz = 0.0, uzz = 0.0;
        tr_uniform_PQ(G_,nu_, 1,0,0, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, uxz,uyz,uzz);

        ASSERT_NEAR(uxy, uyx, tol_);
        ASSERT_NEAR(uyz, -uzy, tol_);
        ASSERT_NEAR(uzx, -uxz, tol_);
    }
}


TEST_F(LiBergerTest, check_triangle_gives_no_NaNs) {
    /*
     * Check the symmetry of the Green's function
     */
     
    double ux = 0.0, uy = 0.0, uz = 0.0;
    tr_uniform_PQ(G_,nu_, 1,1,1, -1,-1, 1,-1, 1,1, 0,0, ux,uy,uz);

    /* Is x is NaN, then x == x is false (or x != x true) */
    ASSERT_TRUE(ux == ux);
    ASSERT_TRUE(uy == uy);
    ASSERT_TRUE(uz == uz);
}


TEST_F(LiBergerTest, triangle_solution_is_rotationally_invariant) {
    /*
     * Create a couple of random triangles and random points. Then rotate
     * triangles and points to see if solution is invariant.
     */
    for (int i = 0; i < n_trials_; i++) {
        double Bx = rand1(), By = rand1();
        double A2x = rand1(), A2y = rand1();
        double A3x = rand1(), A3y = rand1();
        
        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;
        
        tr_uniform_PQ_c(G_, nu_, P, Qx, Qy, Bx, By, A2x, A2y, A3x, A3y,
                        ux, uy, uz);
                        
        for (int i = 0; i < n_trials_; i++) {
            /* Random rotation matrix, this is easy in 2D */
            double phi = rand1()*2*M_PI;
            double s = sin(phi), c = cos(phi);
            
            /* Rotate B, A2, A3 and Q vectors counterclockwise */
            double _Bx = c*Bx - s*By;
            double _By = s*Bx + c*By;
            
            double _A2x = c*A2x - s*A2y;
            double _A2y = s*A2x + c*A2y;
            
            double _A3x = c*A3x - s*A3y;
            double _A3y = s*A3x + c*A3y;
            
            double _Qx = c*Qx - s*Qy;
            double _Qy = s*Qx + c*Qy;

            double _ux = 0.0, _uy = 0.0, _uz = 0.0;
            
            tr_uniform_PQ_c(G_, nu_, P, _Qx, _Qy, _Bx, _By, _A2x, _A2y,
                            _A3x, _A3y, _ux, _uy, _uz);
            
            /* Rotate ux,uy */
            double __ux = c*ux - s*uy;
            double __uy = s*ux + c*uy;
            
            ASSERT_NEAR(_ux, __ux, tol_);
            ASSERT_NEAR(_uy, __uy, tol_);
            ASSERT_NEAR(_uz, uz, tol_);
        }
    }
}


TEST_F(LiBergerTest, order_of_parameters_in_triangle_solution_does_not_matter) {
    /*
     * Check if the order in which the corners of the triangle are specified
     * does matter.
     */
    for (int i = 0; i < n_trials_; i++) {
        double A1x = 2*rand1()-1, A1y = 2*rand1()-1;
        double A2x = 2*rand1()-1, A2y = 2*rand1()-1;
        double A3x = 2*rand1()-1, A3y = 2*rand1()-1;
        double Bx = 10*(2*rand1()-1), By = 10*(2*rand1()-1);
        
        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;

        tr_uniform_PQ(G_,nu_, P,Qx,Qy, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, ux,uy,uz);

        double _ux = 0.0, _uy = 0.0, _uz = 0.0;

        tr_uniform_PQ(G_,nu_, P,Qx,Qy, A2x,A2y, A1x,A1y,
                      A3x,A3y, Bx,By, _ux,_uy,_uz);

        double __ux = 0.0, __uy = 0.0, __uz = 0.0;

        tr_uniform_PQ(G_,nu_, P,Qx,Qy, A2x,A2y, A3x,A3y,
                      A1x,A1y, Bx,By, __ux,__uy,__uz);
                      
        ASSERT_NEAR(_ux, ux, tol_);
        ASSERT_NEAR(_uy, uy, tol_);
        ASSERT_NEAR(_uz, uz, tol_);
        
        ASSERT_NEAR(__ux, ux, tol_);
        ASSERT_NEAR(__uy, uy, tol_);
        ASSERT_NEAR(__uz, uz, tol_);
    }
}


TEST_F(LiBergerTest, triangle_solution_is_scale_invariant) {
    /*
     * Check if we the solution is invariant with respect to rescaling
     */
    for (int i = 0; i < n_trials_; i++) {
        double A1x = 2*rand1()-1, A1y = 2*rand1()-1;
        double A2x = 2*rand1()-1, A2y = 2*rand1()-1;
        double A3x = 2*rand1()-1, A3y = 2*rand1()-1;
        double Bx = 10*(2*rand1()-1), By = 10*(2*rand1()-1);
        
        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;

        tr_uniform_PQ(G_,nu_, P,Qx,Qy, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, ux,uy,uz);
                      
        for (int j = 0; j < n_trials_; j++) {        
            /* scale factor */
            double fac = rand1();
            
            double _ux = 0.0, _uy = 0.0, _uz = 0.0;
            
            tr_uniform_PQ(G_,nu_, P,Qx,Qy, fac*A1x,fac*A1y, fac*A2x,fac*A2y,
                          fac*A3x,fac*A3y, fac*Bx,fac*By, _ux,_uy,_uz);

            ASSERT_NEAR(_ux, fac*ux, tol_);
            ASSERT_NEAR(_uy, fac*uy, tol_);
            ASSERT_NEAR(_uz, fac*uz, tol_);
        }
    }
}


TEST_F(LiBergerTest, triangle_becomes_Boussinesq_Cerruti_at_large_distances) {
    /* We need a lower tolerance here since agreement is only approximate */
    const double tol = 1e-6;

    /*
     * Check if we recover the Boussinesq-Cerruti solution if far away from
     * the homogeneously loaded triangle.
     */
    for (int i = 0; i < n_trials_; i++) {
        double A1x = 2*rand1()-1, A1y = 2*rand1()-1;
        double A2x = 2*rand1()-1, A2y = 2*rand1()-1;
        double A3x = 2*rand1()-1, A3y = 2*rand1()-1;
        double Bx = 100000*(1.0+rand1()), By = 100000*(1.0+rand1());
        
        double area = tr_area(A1x,A1y, A2x,A2y, A3x,A3y);

        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;

        tr_uniform_PQ(G_,nu_, P/area,Qx/area,Qy/area, A1x,A1y, A2x,A2y,
                      A3x,A3y, Bx,By, ux,uy,uz);
                      
        double _ux = 0.0, _uy = 0.0, _uz = 0.0;

        pt_PQ(G_,nu_, P,Qx,Qy, Bx,By,0.0, _ux,_uy,_uz);

        ASSERT_NEAR(_ux, ux, tol);
        ASSERT_NEAR(_uy, uy, tol);
        ASSERT_NEAR(_uz, uz, tol);
    }
}


TEST_F(LiBergerTest, square_solution_from_triangles_agrees_with_analytical) {
    const double tol = 1e-6;
    
    double fac = (1-nu_)/(2*M_PI*G_);

    /*
     * Check if we recover the solution for a uniformly load square
     * (Johnson, p.54, only uz).
     * First check for 1x1 square pushed in the middle, then do random
     * squares.
     */

    double a = 0.5, b = 0.5;
    double x = 0.0, y = 0.0;
        
    double P = 2*rand1()-1, Qx = 0.0, Qy = 0.0;
    
    for (int i = 0; i < n_trials_; i++) {       
        double ux = 0.0, uy = 0.0, uz = 0.0;

        sq_uniform_PQ(G_,nu_, P,Qx,Qy, a,b, x,y, ux,uy,uz);
                      
        double _ux = 0.0, _uy = 0.0, _uz = 0.0;

        _uz = ( (x+a)*log( ( (y+b)+sqrt((y+b)*(y+b)+(x+a)*(x+a)) )/
                           ( (y-b)+sqrt((y-b)*(y-b)+(x+a)*(x+a)) ) )+
                (y+b)*log( ( (x+a)+sqrt((y+b)*(y+b)+(x+a)*(x+a)) ) /
                           ( (x-a)+sqrt((y+b)*(y+b)+(x-a)*(x-a)) ) )+
                (x-a)*log( ( (y-b)+sqrt((y-b)*(y-b)+(x-a)*(x-a)) ) /
                           ( (y+b)+sqrt((y+b)*(y+b)+(x-a)*(x-a)) ) )+
                (y-b)*log( ( (x-a)+sqrt((y-b)*(y-b)+(x-a)*(x-a)) ) /
                           ( (x+a)+sqrt((y-b)*(y-b)+(x+a)*(x+a)) ) ) )*fac*P;

        ASSERT_NEAR(_uz, uz, tol_);
        
        a = rand1(); b = rand1();
        x = 2*rand1()-1; y = 2*rand1()-1;
        
        P = 2*rand1()-1; Qx = 0.0; Qy = 0.0;
    }
}
