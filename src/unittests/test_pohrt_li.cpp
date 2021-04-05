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

#include "boussinesq_cerruti.h"
#include "geometry.h"
#include "li_berger.h"
#include "pohrt_li.h"

#include "gtest/gtest.h"

using namespace PohrtLi;

class PohrtLiTest: public ::testing::Test {
protected:
    PohrtLiTest() {
    };

    virtual ~PohrtLiTest() {
    };
    
    static constexpr int n_trials_ = 100;
    static constexpr double G_ = 1.0;
    static constexpr double nu_ = 1.0/3;

    /* Tolerance for comparing values */
    static constexpr double tol_ = 1e-12;
};


TEST_F(PohrtLiTest, check_square_solution_symmetry) {
    /*
     * Check the symmetry of the Green's function
     */
    for (int i = 0; i < n_trials_; i++) {
	double a = rand1(), b = rand1();
        double Bx = 10*(2*rand1()-1), By = 10*(2*rand1()-1);
       
        double uxx = 0.0, uyx = 0.0, uzx = 0.0;
        sq_uniform_PQ(G_,nu_, 0,1,0, a,b, Bx,By, uxx,uyx,uzx);

        double uxy = 0.0, uyy = 0.0, uzy = 0.0;
        sq_uniform_PQ(G_,nu_, 0,0,1, a,b, Bx,By, uxy,uyy,uzy);

        double uxz = 0.0, uyz = 0.0, uzz = 0.0;
        sq_uniform_PQ(G_,nu_, 1,0,0, a,b, Bx,By, uxz,uyz,uzz);

        ASSERT_NEAR(uxy, uyx, tol_);
        ASSERT_NEAR(uyz, -uzy, tol_);
        ASSERT_NEAR(uzx, -uxz, tol_);
    }
}


TEST_F(PohrtLiTest, square_solution_is_rotationally_invariant) {
    /*
     * Rotate the point of interest by 90, 180, and 270 degrees and see if
     * solution remains invariant.
     */
    for (int i = 0; i < n_trials_; i++) {
        double Bx = 10*rand1(), By = 10*rand1();
        
        //double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        double Qx = 2*rand1()-1, P = 0.0, Qy = 0.0;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;
        
        sq_uniform_PQ(G_, nu_, P, Qx, Qy, 1.0, 1.0, Bx, By, ux, uy, uz);
                        
        for (int i = 0; i < 3; i++) {
            /* Random rotation matrix, this is easy in 2D */
            double phi = (i+1)*2*M_PI/4;
            double s = sin(phi), c = cos(phi);
            
            /* Rotate B, A2, A3 and Q vectors counterclockwise */
            double _Bx = c*Bx - s*By;
            double _By = s*Bx + c*By;
            
            double _Qx = c*Qx - s*Qy;
            double _Qy = s*Qx + c*Qy;

            double _ux = 0.0, _uy = 0.0, _uz = 0.0;
            
            sq_uniform_PQ(G_, nu_, P, _Qx, _Qy, 1.0, 1.0,_Bx, _By,
                          _ux, _uy, _uz);
            
            /* Rotate ux,uy */
            double __ux = c*ux - s*uy;
            double __uy = s*ux + c*uy;
            
            ASSERT_NEAR(_ux, __ux, tol_);
            ASSERT_NEAR(_uy, __uy, tol_);
            ASSERT_NEAR(_uz, uz, tol_);
        }
    }
}


TEST_F(PohrtLiTest, square_solution_is_scale_invariant) {
    /*
     * Check if we the solution is invariant with respect to rescaling
     */
    for (int i = 0; i < n_trials_; i++) {
        double a = 2*rand1()-1, b = 2*rand1()-1;
        double Bx = 10*(2*rand1()-1), By = 10*(2*rand1()-1);
        
        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;

        sq_uniform_PQ(G_,nu_, P,Qx,Qy, a,b, Bx,By, ux,uy,uz);
                      
        for (int j = 0; j < n_trials_; j++) {        
            /* scale factor */
            double fac = rand1();
            
            double _ux = 0.0, _uy = 0.0, _uz = 0.0;
            
            sq_uniform_PQ(G_,nu_, P,Qx,Qy, fac*a,fac*b, fac*Bx,fac*By,
			  _ux,_uy,_uz);

            ASSERT_NEAR(_ux, fac*ux, tol_);
            ASSERT_NEAR(_uy, fac*uy, tol_);
            ASSERT_NEAR(_uz, fac*uz, tol_);
        }
    }
}


TEST_F(PohrtLiTest, square_becomes_Boussinesq_Cerruti_at_large_distances) {
    /* We need a lower tolerance here since agreement is only approximate */
    const double tol = 1e-6;

    /*
     * Check if we recover the Boussinesq-Cerruti solution if far away from
     * the homogeneously loaded triangle.
     */
    for (int i = 0; i < n_trials_; i++) {
        double a = 2*rand1()-1, b = 2*rand1()-1;
        double Bx = 1000000*(1.0+rand1()), By = 1000000*(1.0+rand1());
        
        double area = a*b;

        double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
        
        double ux = 0.0, uy = 0.0, uz = 0.0;

        sq_uniform_PQ(G_,nu_, P/area,Qx/area,Qy/area, a,b, Bx,By, ux,uy,uz);
                      
        double _ux = 0.0, _uy = 0.0, _uz = 0.0;

        pt_PQ(G_,nu_, P,Qx,Qy, Bx,By,0.0, _ux,_uy,_uz);

        ASSERT_NEAR(_ux, ux, tol);
        ASSERT_NEAR(_uy, uy, tol);
        ASSERT_NEAR(_uz, uz, tol);
    }
}


TEST_F(PohrtLiTest, square_solution_equals_Li_Berger) {
    /*
     * Compare to Li-Berger results
     */
    for (double fac = 0.5; fac < 10; fac *= 10) {
        for (int i = 0; i < n_trials_; i++) {
            double Bx = fac*rand1(), By = fac*rand1();

            double P = 2*rand1()-1, Qx = 2*rand1()-1, Qy = 2*rand1()-1;
            double a = rand1(), b = rand1();
        
            double ux = 0.0, uy = 0.0, uz = 0.0;
        
            sq_uniform_PQ(G_,nu_, P,Qx,Qy, a,b, Bx,By, ux,uy,uz);

            double _ux = 0.0, _uy = 0.0, _uz = 0.0;
        
            LiBerger::sq_uniform_PQ(G_,nu_, P,Qx,Qy, a,b, Bx,By,
                                    _ux,_uy,_uz); 

            ASSERT_NEAR(ux, _ux, tol_);
            ASSERT_NEAR(uy, _uy, tol_);
            ASSERT_NEAR(uz, _uz, tol_);
        }
    }
}
