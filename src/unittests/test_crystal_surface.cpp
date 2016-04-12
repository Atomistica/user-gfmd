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

#include "error.h"
#include "fcc100_surface.h"

#include "gtest/gtest.h"

class CrystalSurfaceTest: public ::testing::Test {
protected:
    CrystalSurfaceTest() {
    };

    virtual ~CrystalSurfaceTest() {
    };

	Error error_;
};


/*!
 * Convert distance to neighbor shell number
 */
int dist_to_neighbor_shell(double dist2) {
    double bordermax = 2.0 /1.1; 
    double border2to3 = 1.6 /1.1;
    double border1to2 = 1.3 /1.1;
    double bordermin = 1.0 /1.1;

    if (dist2 > bordermax * bordermax) return -1;
    else if (dist2 > border2to3 * border2to3) return 2;
    else if (dist2 > border1to2 * border1to2) return 1;
    else if (dist2 > bordermin * bordermin) return 0;
    return -1;
}


TEST_F(CrystalSurfaceTest, fcc100) {
	FCC100Surface surface(1.0, 1, &error_);

    // Construct neighbor list up to third neighbors.
    surface.compute_neighbor_list(1.8);

    // Check if all neighbors actually point to atoms.
	ASSERT_TRUE(surface.check_neighbors(&error_));
}


TEST_F(CrystalSurfaceTest, fcc100_supercell) {
    for (int i = 0; i < 2; i++) {
    	FCC100Surface surface(1.0, 1, &error_);

        // Construct neighbor list up to third neighbors.
        surface.compute_neighbor_list(1.8);

        // Check if all neighbors actually point to atoms.
	    ASSERT_TRUE(surface.check_neighbors(&error_));

        // Copy surface.
        FCC100Surface surface2(surface);
        FCC100Surface surface3(surface);

        // Create supercell.
        if (i == 0) {
            surface.supercell(2, 2, 2, &error_);
        }
        else {
            surface.supercell(2, 2, 2, &error_, 0, 1, 0);
        };

        // Check if all neighbors actually point to atoms.
	   ASSERT_TRUE(surface.check_neighbors(&error_));

        // Create supercell of second surface, and create virgin neighbor list.
        if (i == 0) {
            surface2.supercell(2, 2, 2, &error_);
        }
        else {
            surface2.supercell(2, 2, 2, &error_, 0, 1, 0);
        };
        surface2.compute_neighbor_list(1.8);

        // Check if all neighbors actually point to atoms.
    	ASSERT_TRUE(surface2.check_neighbors(&error_));

        // Check if neighbors list of surface and surface2 agree.
        const double tol = 1e-6;
        ASSERT_EQ(surface.get_number_of_neighbors(), surface2.get_number_of_neighbors());
        const CrystalSurface::lattice_neighbor_t *neigh1 = surface.get_neighbors();
        for (int i = 0; i < surface.get_number_of_neighbors(); i++) {
            const CrystalSurface::lattice_neighbor_t *neigh2 = surface2.get_neighbors();
            for (int j = 0; j< surface2.get_number_of_neighbors(); j++) {
                if (neigh1->indexi == neigh2->indexi &&
                    neigh1->indexj == neigh2->indexj &&
                    neigh1->x == neigh2->x &&
                    neigh1->y == neigh2->y &&
                    neigh1->ab == neigh2->ab) {
                    ASSERT_TRUE(std::abs(neigh1->rnorm-neigh2->rnorm) < tol);
                    ASSERT_TRUE(std::abs(neigh1->r[0]-neigh2->r[0]) < tol);
                    ASSERT_TRUE(std::abs(neigh1->r[1]-neigh2->r[1]) < tol);
                    ASSERT_TRUE(std::abs(neigh1->r[2]-neigh2->r[2]) < tol);
                }
                neigh2++;
            }

            neigh1++;
        }
    }
}


TEST_F(CrystalSurfaceTest, reverse_indices) {
	FCC100Surface surface(1.0, 1, &error_);

    // Construct neighbor list up to third neighbors.
    surface.compute_neighbor_list(1.8);

    // Check if all neighbors actually point to atoms.
	ASSERT_TRUE(surface.check_neighbors(&error_));

    // Reverse indices.
    surface.reverse_indices();

    // Check if all neighbors actually point to atoms.
	ASSERT_TRUE(surface.check_neighbors(&error_));

    // Create supercell.
    surface.supercell(2, 2, 2, &error_);

    // Check all neighbors
    ASSERT_TRUE(surface.check_neighbors(&error_));
}


TEST_F(CrystalSurfaceTest, fcc100_3nn) {
    // Check all neighbors up to third

    // Distance between alpha=n and alpha=n+1 (single atom unit cells)
    double sqrt2 = sqrt(2.0);
    double invsqrt2 = 1.0/sqrt2;

    const int NUMLNS = 42;
    double ref_r[3*NUMLNS] = {  
        // ab = 0
         1, 0,0, 
         0, 1,0,
        -1, 0,0,
         0,-1,0,
         1, 1,0,
         1,-1,0,
        -1, 1,0,
        -1,-1,0,
        // ab = 1
         0.5, 0.5, invsqrt2, 
        -0.5, 0.5, invsqrt2,
         0.5,-0.5, invsqrt2,
        -0.5,-0.5, invsqrt2,
         1.5, 0.5, invsqrt2,
         0.5, 1.5, invsqrt2,
        -1.5, 0.5, invsqrt2,
        -0.5, 1.5, invsqrt2,
         1.5,-0.5, invsqrt2,
         0.5,-1.5, invsqrt2,
        -1.5,-0.5, invsqrt2,
        -0.5,-1.5, invsqrt2,
        // ab = -1
         0.5, 0.5,-invsqrt2,
        -0.5, 0.5,-invsqrt2,
         0.5,-0.5,-invsqrt2,
        -0.5,-0.5,-invsqrt2,
         1.5, 0.5,-invsqrt2,
         0.5, 1.5,-invsqrt2,
        -1.5, 0.5,-invsqrt2,
        -0.5, 1.5,-invsqrt2,
         1.5,-0.5,-invsqrt2,
         0.5,-1.5,-invsqrt2,
        -1.5,-0.5,-invsqrt2,
        -0.5,-1.5,-invsqrt2,
        // ab = 2 
         0, 0, sqrt2,
         1, 0, sqrt2,
         0, 1, sqrt2,
        -1, 0, sqrt2,
         0,-1, sqrt2,
        // ab = -2
         0, 0,-sqrt2,
         1, 0,-sqrt2,
         0, 1,-sqrt2,
        -1, 0,-sqrt2,
         0,-1,-sqrt2
    };
    bool found[NUMLNS];

    std::fill(found, found+NUMLNS, false);

	FCC100Surface surface(1.0, 1, &error_);

    // Construct neighbor list up to third neighbors.
    surface.compute_neighbor_list(1.8);

    const CrystalSurface::lattice_neighbor_t *neigh = surface.get_neighbors(0);
    for (int i = 0; i < surface.get_number_of_neighbors(0); i++) {
        for (int j = 0; j < NUMLNS; j++) {
            vec<double,3> r(neigh->r);
            r -= &ref_r[3*j];
            if (r.nrm2() < 1e-6) {
                ASSERT_FALSE(found[j]);
                found[j] = true;
            }
        }
        neigh++;
    }

    for (int i = 0; i < NUMLNS; i++) {
        ASSERT_TRUE(found[i]);
    }
}


TEST_F(CrystalSurfaceTest, fcc100_neighbor_shells) {
	FCC100Surface surface(1.0, 1, &error_);

    // Construct neighbor list up to third neighbors.
    surface.compute_neighbor_list(1.8);

    // Check if neighbor shells are okay.
    const CrystalSurface::lattice_neighbor_t *neigh = surface.get_neighbors(0);
    for (int i = 0; i < surface.get_number_of_neighbors(0); i++) {
        ASSERT_TRUE(neigh->neighbor_shell == 
                    dist_to_neighbor_shell(neigh->rnorm*neigh->rnorm));
        neigh++;
    }
}
