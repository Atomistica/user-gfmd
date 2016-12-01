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
#ifndef __NONPERIODIC_STIFFNESS_H
#define __NONPERIODIC_STIFFNESS_H

#include <fftw3.h>

#include "domain.h"
#include "error.h"
#include "linearalgebra.h"
#include "memory.h"
#include "surface_stiffness.h"

namespace LAMMPS_NS {

class NonperiodicStiffnessKernel : public StiffnessKernel {
 public:
    NonperiodicStiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
    virtual ~NonperiodicStiffnessKernel();
  
    virtual void get_stiffness_matrix(double, double, double_complex *,
                                      double_complex dU);
    virtual void get_stiffness_matrix(int, int, double, double, double_complex *,
                                      double_complex dU); 

 private:
    /*!
     * Crystal lattice
     */
    CrystalSurface *crystal_surface_;

    /*!
     * Poisson number
     */
    double Poisson_number_;

    /*!
     * Shear modulus
     */
    double shear_modulus_;

    /*!
     * Grid dimension
     */
    int nx_, ny_;

    /*!
     * FFT buffer
     */
    fftw_complex **complex_buffer_;

    void create_and_fill_buffer(int, int);
};

}

#endif
