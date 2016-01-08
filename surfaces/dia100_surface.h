#ifdef CRYSTAL_SURFACE_CLASS

CrystalSurfaceStyle(dia100,Diamond100Surface)

#else

#ifndef __DIAMOND100_SURFACE_H
#define __DIAMOND100_SURFACE_H

#include "crystal_surface.h"

class Diamond100Surface : public CrystalSurface {
 public:
	Diamond100Surface(int, int *, char **, Error *);
	Diamond100Surface(double, int, Error *);

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
