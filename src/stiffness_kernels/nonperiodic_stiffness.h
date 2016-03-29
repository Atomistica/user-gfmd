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
    ~NonperiodicStiffnessKernel();
  
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
