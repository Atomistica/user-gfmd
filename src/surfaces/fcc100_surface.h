#ifdef CRYSTAL_SURFACE_CLASS

CrystalSurfaceStyle(fcc100,FCC100Surface)

#else

#ifndef __FCC100_SURFACE_H
#define __FCC100_SURFACE_H

#include "crystal_surface.h"

class FCC100Surface : public CrystalSurface {
 public:
	FCC100Surface(int, int *, char **, Error *);
	FCC100Surface(double, int, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *f);

 protected:
    /*!
     * Create primitive unit cell
     */
    void primitive_cell(Error *error);

    /*!
     * Number of unit cells (along surface normal).
     */
    int nu_;

    /*!
     * The lattice constant.
     */
    double lattice_constant_;
};

#endif

#endif
