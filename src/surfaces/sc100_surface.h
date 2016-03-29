#ifdef CRYSTAL_SURFACE_CLASS

CrystalSurfaceStyle(sc100,SC100Surface)

#else

#ifndef __SC100_SURFACE_H
#define __SC100_SURFACE_H

#include "crystal_surface.h"

class SC100Surface : public CrystalSurface {
 public:
	SC100Surface(int, int *, char **, Error *);
	SC100Surface(double, int, Error *);

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
